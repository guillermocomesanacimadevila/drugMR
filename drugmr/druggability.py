import time
from datetime import date

import polars as pl
import pymysql
import requests
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings

# the ChEMBL client fetches 20 rows per request by default; EGFR activity alone is ~17k rows.
# 500 per page was faster than 1000 in a test (85 s vs 133 s for EGFR + ACE)
Settings.Instance().MAX_LIMIT = 500


class DruggableGenome:

    """ 
    Query targets to Pharos database 
    Tclin enrichment -> background.tsv plus a Fisher test (scipy)
    """

    def __init__(
            self,
            database="pharos319",
            host="tcrd.ncats.io",
            user="tcrd"
    ):

        self.database = database
        self.connection = pymysql.connect(
            host=host,
            user=user,
            password="",
            database=database
        )

    def query(self):
        cursor = self.connection.cursor()

        sql = """
SELECT p.uniprot,
       p.sym,
       x.value AS ensembl_gene_id,
       t.tdl,
       t.fam
FROM target t
JOIN t2tc ON t2tc.target_id = t.id
JOIN protein p ON p.id = t2tc.protein_id
LEFT JOIN xref x ON x.protein_id = p.id
                AND x.xtype = 'Ensembl'
                AND x.value LIKE 'ENSG%';
        """
        cursor.execute(sql)
        rows = cursor.fetchall()
        
        ###### check
        # print(rows)

        columns = []
        for column in cursor.description:
            columns.append(column[0])

        cursor.close()
        return pl.DataFrame(rows, schema=columns, orient="row")

    def close(self):
        self.connection.close()

    def find_duplicates(self, table):
        counts = table.group_by("ensembl_gene_id").agg(pl.len().alias("n_rows"))
        duplicated = counts.filter(pl.col("n_rows") > 1)
        return duplicated

    ##### --- --- --- --- --- --- --- --- ---
    ###### Focus on functions below this! ***
    ##### --- --- --- --- --- --- --- --- ---
    
    def collapse_duplicates(self, table):
        tdl_rank = {"Tclin": 1, "Tchem": 2, "Tbio": 3, "Tdark": 4}

        # rank each row so the most developed TDL sorts first
        table = table.with_columns(pl.col("tdl").replace_strict(tdl_rank, default=5).alias("tdl_rank"))
        table = table.sort(["ensembl_gene_id", "tdl_rank"])

        # keep the first (best) row per Ensembl ID
        collapsed = table.group_by("ensembl_gene_id", maintain_order=True).first()

        # flag Ensembl IDs that matched more than one Pharos row
        duplicated = self.find_duplicates(table)
        collapsed = collapsed.with_columns(
            pl.col("ensembl_gene_id")
            .is_in(duplicated["ensembl_gene_id"].implode())
            .alias("pharos_multiple_matches")
        )

        collapsed = collapsed.drop("tdl_rank")
        return collapsed

    def run(self, ensembl_ids):
        pharos = self.query()
        pharos = pharos.filter(
            pl.col("ensembl_gene_id").is_not_null()
        )
        pharos = self.collapse_duplicates(pharos)
        targets = pl.DataFrame(
            {
                "ensembl_gene_id": ensembl_ids
            }
        )
        result = targets.join(pharos, on="ensembl_gene_id", how="left")
        result = result.with_columns(
            pl.col("tdl").is_null().alias("not_in_pharos"),
            pl.lit(self.database).alias("pharos_release"),
        )
        return result


class KnownDrugs:
    
    """ 
    CheMBL REST API
    - Gather:
    -> Known drugs to target X
    -> Direction of drug to target X
    -> Disease which drug Y was designed for 

    # KnownDrugs brain access -> CNS MPO score
    # Potent compounds (pChEMBL ≥ 6)
    """

    # ChEMBL action_type -> what the drug does to the target protein.
    # Covers every action_type in ChEMBL_37 mechanisms; anything else is "Unclassified"
    # (BINDING AGENT, MODULATOR, VACCINE ANTIGEN, enzymes, chelating/cross-linking agents, OTHER...)
    action_direction = {
        "INHIBITOR": "Inhibit",
        "ANTAGONIST": "Inhibit",
        "ALLOSTERIC ANTAGONIST": "Inhibit",
        "BLOCKER": "Inhibit",
        "NEGATIVE ALLOSTERIC MODULATOR": "Inhibit",
        "NEGATIVE MODULATOR": "Inhibit",
        "INVERSE AGONIST": "Inhibit",
        "ANTISENSE INHIBITOR": "Inhibit",
        "RNAI INHIBITOR": "Inhibit",
        "GENE EDITING NEGATIVE MODULATOR": "Inhibit",
        "DEGRADER": "Inhibit",
        "DISRUPTING AGENT": "Inhibit",
        "SEQUESTERING AGENT": "Inhibit",
        "AGONIST": "Activate",
        "PARTIAL AGONIST": "Activate",
        "POSITIVE ALLOSTERIC MODULATOR": "Activate",
        "POSITIVE MODULATOR": "Activate",
        "ACTIVATOR": "Activate",
        "OPENER": "Activate",
        "STABILISER": "Activate",
        "EXOGENOUS PROTEIN": "Activate",
        "EXOGENOUS GENE": "Activate",
    }

    def __init__(self):
        self.target = new_client.target
        self.molecule = new_client.molecule
        self.mechanism = new_client.mechanism
        self.indication = new_client.drug_indication
        self.activity = new_client.activity
        self.release = self.get_release()

        # B3DB: blood brain barrier labels, pinned to one commit of github.com/theochem/B3DB
        self.b3db_release = "B3DB 75dab1cc (2025-09-10)"
        self.b3db_url = "https://raw.githubusercontent.com/theochem/B3DB/75dab1cc607bb7a03f3de5c576cffd223e767844/B3DB/B3DB_classification.tsv"
        b3db = pl.read_csv(self.b3db_url, separator="\t", infer_schema_length=None)
        self.b3db = b3db.select([
            pl.col("Inchi").alias("standard_inchi"), # same name as the ChEMBL column, for the join
            pl.col("compound_name").alias("b3db_name"),
            pl.col("BBB+/BBB-").alias("b3db_label"),
            pl.col("logBB").alias("b3db_logbb"),
        ])

    def get_release(self):
        url = "https://www.ebi.ac.uk/chembl/api/data/status.json"
        for attempt in range(3):
            response = requests.get(url)
            if response.status_code == 200:
                return response.json()["chembl_db_version"]
            print(f"ChEMBL status check failed (try {attempt + 1} of 3), retrying")
            time.sleep(5)
        response.raise_for_status()

    @staticmethod
    def chembl_to_df(client_self):
        return pl.DataFrame(list(client_self), infer_schema_length=None)

    def query_target_table(
            self,
            uniprot_targets,
            organism: str = "Homo sapiens"
        ):
        # single proteins plus groups that contain them (families, complexes),
        # because some drugs are recorded against the group, e.g. pan-HER inhibitors
        target_types = [
            "SINGLE PROTEIN",
            "PROTEIN FAMILY",
            "PROTEIN COMPLEX",
            "PROTEIN COMPLEX GROUP",
            "SELECTIVITY GROUP",
        ]
        query = self.target.filter(
            target_components__accession__in=uniprot_targets,
            target_type__in=target_types,
            organism=organism,
        ).only(["target_chembl_id", "target_type", "target_components"]) # only these cols
        df = self.chembl_to_df(query)

        # none of these proteins are in ChEMBL: return an empty table with the right columns
        if df.is_empty():
            return pl.DataFrame(schema={
                "target_chembl_id": pl.String,
                "accession": pl.String,
                "target_type": pl.String,
                "via_group": pl.Boolean,
            })

        df = df.explode("target_components", empty_as_null=True)
        df = df.unnest("target_components")

        # a family or complex lists all its member proteins; keep only the rows for our targets
        df = df.filter(pl.col("accession").is_in(uniprot_targets))

        # flag drugs that will come through a group rather than the protein alone
        df = df.with_columns(
            (pl.col("target_type") != "SINGLE PROTEIN").alias("via_group")
        )
        df = df.drop(["component_id", "target_component_synonyms", "target_component_xrefs"])
        return df

    def query_molecule_table(
            self,
            chembl_ids,
    ):
        query = self.molecule.filter(
            molecule_chembl_id__in=chembl_ids, # molecule_chembl_id col
        )
        df = self.chembl_to_df(query)

        # no drugs to look up: return an empty table with the right columns
        if df.is_empty():
            return pl.DataFrame(schema={
                "molecule_chembl_id": pl.String,
                "pref_name": pl.String,
                "molecule_type": pl.String,
                "first_approval": pl.Int64,
                "withdrawn_flag": pl.Boolean,
                "black_box_warning": pl.Int64,
                "chemical_probe": pl.Int64,
                "oral": pl.Boolean,
                "parenteral": pl.Boolean,
                "alogp": pl.Float64,
                "mw_freebase": pl.Float64,
                "psa": pl.Float64,
                "hbd": pl.Int64,
                "num_ro5_violations": pl.Int64,
                "standard_inchi_key": pl.String,
                "standard_inchi": pl.String,
                "canonical_smiles": pl.String,
            })

        # biologics have no properties or structures: if every molecule is one, polars reads these as Null, not Struct
        property_columns = {
            "alogp": pl.Float64,
            "mw_freebase": pl.Float64,
            "psa": pl.Float64,
            "hbd": pl.Int64,
            "num_ro5_violations": pl.Int64,
        }
        structure_columns = {
            "standard_inchi_key": pl.String,
            "standard_inchi": pl.String,
            "canonical_smiles": pl.String,
        }
        if df.schema["molecule_properties"] == pl.Null:
            df = df.drop("molecule_properties")
            for column, dtype in property_columns.items():
                df = df.with_columns(pl.lit(None, dtype=dtype).alias(column))
        else:
            df = df.unnest("molecule_properties")
        if df.schema["molecule_structures"] == pl.Null:
            df = df.drop("molecule_structures")
            for column, dtype in structure_columns.items():
                df = df.with_columns(pl.lit(None, dtype=dtype).alias(column))
        else:
            df = df.unnest("molecule_structures")
        df = df.select([
            "molecule_chembl_id",
            "pref_name",
            "molecule_type",
            "first_approval",
            "withdrawn_flag",
            "black_box_warning", # binary -> 1 == Y - 0 == N
            # "therapeutic_flag"
            "chemical_probe", # binary -> #### if 0 -> compound == drug -> if 1 then compound is a chemical probe
            # chemical probe -> research-based -> used in cells / animals -> knock protein -> optimised for potency and selectivity -> just to show an effect in lab
            # chemical probe == 0 -> drug (used in humans != in ["cells", "animals"])
            "oral",
            ##### -> logD/pKa
            "parenteral",
            pl.col("alogp").cast(pl.Float64), # lipophilicity

            # Smaller molecules cross membranes and the blood brain barrier more easily
            # CNS drugs are usually under about 450 Da, and in CNS MPO the score is best at ≤ 360
            pl.col("mw_freebase").cast(pl.Float64), # molecular weight (in daltons) -> if this is low + CNS MPO == high then high chance of it working for CNS!!!!
            
            # BBB reach under 60-70 angstroms -> under 120-140 (good oral permeability and bioabailability)
            # if greater than 140 struggles to surpass lipid cell membranes 
            pl.col("psa").cast(pl.Float64), # sum of the surface area of all polar atoms within a molecule
            "hbd", # hydrogen bond donor
            "num_ro5_violations", # how many of Lipinski's 4 oral drug rules it breaks (0 to 4); 0 ideal, >= 2 poor oral absorption
            "standard_inchi_key", # structure-related annotation
            "standard_inchi", # full InChI, used to match B3DB
            "canonical_smiles", # structure-related annotation
        ])
        return df

    def query_mechanism_table(
            self,
            target_chembl_ids, # the moelcules itself -> left join with target_chembl_ids()
    ):
        query = self.mechanism.filter(
            target_chembl_id__in=target_chembl_ids,
        )
        df = self.chembl_to_df(query)

        # no drugs for any of these targets: return an empty table with the right columns
        if df.is_empty():
            return pl.DataFrame(schema={
                "molecule_chembl_id": pl.String,
                "parent_molecule_chembl_id": pl.String,
                "target_chembl_id": pl.String,
                "mechanism_of_action": pl.String,
                "action_type": pl.String,
                "max_phase": pl.Int64,
            })

        df = df.select([
            "molecule_chembl_id",
            "parent_molecule_chembl_id",
            "target_chembl_id",
            "mechanism_of_action",
            "action_type",
            "max_phase",
        ])

        # salts and the parent drug are separate molecules; collapse them onto the parent
        # and keep the highest phase per target, parent drug and action type
        df = df.sort("max_phase", descending=True)
        df = df.unique(
            subset=["target_chembl_id", "parent_molecule_chembl_id", "action_type"],
            keep="first",
            maintain_order=True,
        )
        return df

    def query_indication_table(
            self,
            drug_ids,

    ):

        query = self.indication.filter(
            parent_molecule_chembl_id__in=drug_ids,
        )
        df = self.chembl_to_df(query)

        # no indications for any of these drugs: return an empty table with the right columns
        if df.is_empty():
            return pl.DataFrame(schema={
                "molecule_chembl_id": pl.String,
                "parent_molecule_chembl_id": pl.String,
                "efo_id": pl.String,
                "efo_term": pl.String,
                "max_phase_for_ind": pl.Float64,
            })

        # print(df.columns)
        df = df.select([
            "molecule_chembl_id",
            "parent_molecule_chembl_id",
            "efo_id",
            "efo_term",
            pl.col("max_phase_for_ind").cast(pl.Float64), # comes back as a string, e.g. "4.0"
        ])
        return df

    def query_activity_table(
            self,
            target_chembl_ids,
            pchembl_thresh,
    ):
        query = self.activity.filter(
            target_chembl_id__in=target_chembl_ids,
            pchembl_value__gte=pchembl_thresh,
        ).only([
            "target_chembl_id",
            "molecule_chembl_id",
            "standard_type",
            "pchembl_value",
        ]) # only these cols
        df = self.chembl_to_df(query)

        # no potent compounds for any of these targets: return an empty table with the right columns
        if df.is_empty():
            return pl.DataFrame(schema={
                "target_chembl_id": pl.String,
                "molecule_chembl_id": pl.String,
                "standard_type": pl.String,
                "pchembl_value": pl.Float64,
            })

        df = df.select([
            "target_chembl_id",
            "molecule_chembl_id",
            "standard_type",
            pl.col("pchembl_value").cast(pl.Float64),
        ])
        return df

    def run(
            self,
            target_list: list[str], # asuming uniprot ids
            pchembl_thresh=6,
        ):

        """ Function in charge of runnin the entire CheMBL-based drugagbility checks """

        # targets table 
        targets = self.query_target_table(target_list)

        # mechanism table
        target_ids = targets["target_chembl_id"].unique().to_list()
        mechanisms = self.query_mechanism_table(target_ids)

        # activities table: single proteins only, a family or selectivity group says little
        # about whether a compound hits this particular protein
        single_ids = targets.filter(~pl.col("via_group"))["target_chembl_id"].unique().to_list()
        activities = self.query_activity_table(single_ids, pchembl_thresh=pchembl_thresh)

        # molecules table
        drug_ids = mechanisms["parent_molecule_chembl_id"].unique().to_list()
        molecules = self.query_molecule_table(drug_ids)

        # add blood brain barrier labels; drugs not in B3DB (and antibodies) get nulls
        molecules = molecules.join(self.b3db, on="standard_inchi", how="left")

        # indications table
        indications = self.query_indication_table(drug_ids)

        return {
            "targets": targets,
            "mechanisms": mechanisms,
            "activities": activities,
            "molecules": molecules,
            "indications": indications,
        } # Brain access for drugs (B3DB)

    def summarise_targets(self, chembl):
        """ One row per UniProt accession: drug counts, phase, directions, potent compounds """

        targets = chembl["targets"].select(["target_chembl_id", "accession", "via_group"])

        # drugs: each mechanism row gets its protein (accession) and direction
        drugs = chembl["mechanisms"].join(targets, on="target_chembl_id", how="inner")
        drugs = drugs.with_columns(
            pl.col("action_type")
            .replace_strict(self.action_direction, default="Unclassified")
            .alias("drug_direction")
        )
        drug_summary = drugs.group_by("accession").agg(
            pl.col("parent_molecule_chembl_id").n_unique().alias("n_drugs"),
            pl.col("parent_molecule_chembl_id").filter(pl.col("max_phase") == 4).n_unique().alias("n_drugs_approved"),
            pl.col("max_phase").max().alias("max_phase"),
            pl.col("parent_molecule_chembl_id").filter(pl.col("drug_direction") == "Inhibit").n_unique().alias("n_drugs_inhibit"),
            pl.col("parent_molecule_chembl_id").filter(pl.col("drug_direction") == "Activate").n_unique().alias("n_drugs_activate"),
            pl.col("via_group").all().alias("drugs_via_group_only"),
        )

        # potent compounds (activity is already single proteins only)
        activities = chembl["activities"].join(targets, on="target_chembl_id", how="inner")
        activity_summary = activities.group_by("accession").agg(
            pl.col("molecule_chembl_id").n_unique().alias("n_potent_compounds"),
            pl.col("pchembl_value").max().alias("best_pchembl"),
        )

        # every accession ChEMBL knows, with or without drugs or potent compounds
        summary = targets.select("accession").unique()
        summary = summary.join(drug_summary, on="accession", how="left")
        summary = summary.join(activity_summary, on="accession", how="left")
        return summary

    def drug_table(self, chembl):
        """ One row per UniProt accession x parent drug: direction, phase, safety flags, brain access, indications """

        targets = chembl["targets"].select(["target_chembl_id", "accession", "via_group"])

        drugs = chembl["mechanisms"].join(targets, on="target_chembl_id", how="inner")
        drugs = drugs.with_columns(
            pl.col("action_type")
            .replace_strict(self.action_direction, default="Unclassified")
            .alias("drug_direction")
        )

        # a drug can reach the same protein through several ChEMBL targets (single + family)
        # or with two action types; collapse to one row per protein and drug
        drugs = drugs.group_by(["accession", "parent_molecule_chembl_id"]).agg(
            pl.col("action_type").unique().sort().str.join("; ").alias("action_type"),
            pl.col("drug_direction").unique().sort().alias("directions"),
            pl.col("max_phase").max().alias("max_phase"),
            pl.col("via_group").all().alias("via_group"),
            pl.col("mechanism_of_action").first().alias("mechanism_of_action"),
            pl.col("target_chembl_id").unique().sort().str.join("; ").alias("target_chembl_id"),
        )

        # one direction per drug: Mixed if it is recorded as both, Unclassified only if nothing else
        drugs = drugs.with_columns(
            pl.when(pl.col("directions").list.contains("Inhibit") & pl.col("directions").list.contains("Activate"))
            .then(pl.lit("Mixed"))
            .when(pl.col("directions").list.contains("Inhibit"))
            .then(pl.lit("Inhibit"))
            .when(pl.col("directions").list.contains("Activate"))
            .then(pl.lit("Activate"))
            .otherwise(pl.lit("Unclassified"))
            .alias("drug_direction")
        )
        drugs = drugs.drop("directions")

        # drug properties and B3DB label
        molecules = chembl["molecules"].select([
            pl.col("molecule_chembl_id").alias("parent_molecule_chembl_id"),
            pl.col("pref_name").alias("drug_name"),
            "molecule_type",
            "first_approval",
            pl.col("withdrawn_flag").alias("withdrawn"),
            (pl.col("black_box_warning") == 1).alias("black_box"),
            "oral",
            "parenteral",
            "alogp",
            "mw_freebase",
            "psa",
            "hbd",
            "num_ro5_violations",
            "standard_inchi_key",
            "canonical_smiles",
            "b3db_label",
            "b3db_logbb",
        ])
        drugs = drugs.join(molecules, on="parent_molecule_chembl_id", how="left")

        # brain access: B3DB where it has the drug; not applicable for biologics
        # (CNS MPO fallback for small molecules not in B3DB still to be decided)
        drugs = drugs.with_columns(
            pl.when(pl.col("molecule_type") != "Small molecule")
            .then(pl.lit("not applicable"))
            .when(pl.col("b3db_label").is_not_null())
            .then(pl.col("b3db_label"))
            .otherwise(pl.lit("unknown"))
            .alias("brain_access")
        )

        # indications: approved ones as text, plus an Alzheimer flag (any phase)
        indications = chembl["indications"]
        approved = indications.filter(pl.col("max_phase_for_ind") == 4).group_by("parent_molecule_chembl_id").agg(
            pl.col("efo_term").unique().sort().str.join("; ").alias("approved_indications")
        )
        alzheimer = indications.group_by("parent_molecule_chembl_id").agg(
            pl.col("efo_term").str.to_lowercase().str.contains("alzheimer").any().alias("ad_indication")
        )
        # every indication with the highest phase reached for it, highest phase first
        all_indications = indications.drop_nulls("efo_term").group_by(["parent_molecule_chembl_id", "efo_term"]).agg(
            pl.col("max_phase_for_ind").max().alias("phase")
        )
        all_indications = all_indications.sort(["phase", "efo_term"], descending=[True, False])
        all_indications = all_indications.with_columns(
            pl.when(pl.col("phase") == 0.5)
            .then(pl.lit("0.5"))
            .otherwise(pl.col("phase").cast(pl.Int64).cast(pl.String))
            .alias("phase")
        )
        all_indications = all_indications.with_columns(
            (pl.col("efo_term") + " (" + pl.col("phase") + ")").alias("indication")
        )
        all_indications = all_indications.group_by("parent_molecule_chembl_id", maintain_order=True).agg(
            pl.col("indication").str.join("; ").alias("indications")
        )
        drugs = drugs.join(approved, on="parent_molecule_chembl_id", how="left")
        drugs = drugs.join(all_indications, on="parent_molecule_chembl_id", how="left")
        drugs = drugs.join(alzheimer, on="parent_molecule_chembl_id", how="left")
        drugs = drugs.with_columns(pl.col("ad_indication").fill_null(False))

        drugs = drugs.rename({"parent_molecule_chembl_id": "drug_chembl_id"})
        drugs = drugs.sort(["accession", "max_phase", "drug_name"], descending=[False, True, False], nulls_last=True)
        return drugs
        

class Tractability: 

    """
    OpenTargets API
    ## subcellularLocations, targetClass, baselineExpression
    - https://api.platform.opentargets.org/api/v4/graphql
    - base_url = "https://api.platform.opentargets.org/api/v4/graphql"
    """

    def __init__(
            self,
            url="https://api.platform.opentargets.org/api/v4/graphql",
            batch_size=25,
    ):
        self.url = url
        self.batch_size = batch_size # Ensembl IDs per request
        self.release = self.get_release()

    def get_release(self):
        query = "{ meta { dataVersion { year month } } }"
        for attempt in range(3):
            response = requests.post(self.url, json={"query": query})
            if response.status_code == 200:
                version = response.json()["data"]["meta"]["dataVersion"]
                return version["year"] + "." + version["month"] # e.g. "26.09"
            print(f"Open Targets release check failed (try {attempt + 1} of 3), retrying")
            time.sleep(5)
        response.raise_for_status()

    def fetch(
            self,
            target_list: list[str], # Ensembl gene IDs
        ):

        graphql_query = """
        query tractability($ids: [String!]!) {
          targets(ensemblIds: $ids) {
            id
            approvedSymbol
            tractability {
              modality
              label
              value
            }
          }
        }
        """

        rows = []

        # send the Ensembl IDs in batches of self.batch_size
        for start in range(0, len(target_list), self.batch_size):
            batch = target_list[start:start + self.batch_size]
            response = requests.post(
                self.url,
                json={"query": graphql_query, "variables": {"ids": batch}},
            )
            response.raise_for_status()
            data = response.json()
            if "errors" in data:
                raise RuntimeError(data["errors"])

            # one row per target x modality x label
            for target in data["data"]["targets"]:
                for item in target["tractability"]:
                    rows.append({
                        "ensembl_gene_id": target["id"],
                        "symbol": target["approvedSymbol"],
                        "modality": item["modality"],
                        "label": item["label"],
                        "value": item["value"],
                    })

        df = pl.DataFrame(rows, schema={
            "ensembl_gene_id": pl.String,
            "symbol": pl.String,
            "modality": pl.String,
            "label": pl.String,
            "value": pl.Boolean,
        })
        return df

    def summarise(self, tractability):
        """ One row per Ensembl ID: True if any flag in the group is True (PROTAC left out) """

        clinical = ["Approved Drug", "Advanced Clinical", "Phase 1 Clinical"]
        groups = [
            ("ot_sm_clinical", "SM", clinical),
            ("ot_sm_discovery", "SM", ["Structure with Ligand", "High-Quality Ligand", "High-Quality Pocket", "Med-Quality Pocket", "Druggable Family"]),
            ("ot_ab_clinical", "AB", clinical),
            ("ot_ab_location", "AB", ["UniProt loc high conf", "GO CC high conf", "Human Protein Atlas loc"]),
            ("ot_oc_clinical", "OC", clinical),
        ]

        summary = tractability.select("ensembl_gene_id").unique()
        for name, modality, labels in groups:
            part = tractability.filter(
                (pl.col("modality") == modality) & pl.col("label").is_in(labels)
            )
            part = part.group_by("ensembl_gene_id").agg(pl.col("value").any().alias(name))
            summary = summary.join(part, on="ensembl_gene_id", how="left")

        # every True label per modality, strongest first (Open Targets order)
        label_order = {
            "SM": clinical + ["Structure with Ligand", "High-Quality Ligand", "High-Quality Pocket", "Med-Quality Pocket", "Druggable Family"],
            "AB": clinical + ["UniProt loc high conf", "GO CC high conf", "UniProt loc med conf", "UniProt SigP or TMHMM", "GO CC med conf", "Human Protein Atlas loc"],
            "OC": clinical,
        }
        label_columns = {
            "SM": "ot_small_molecule_labels",
            "AB": "ot_antibody_labels",
            "OC": "ot_oligonucleotide_labels",
        }
        sm_tractable = clinical + ["Structure with Ligand", "High-Quality Ligand"]
        ab_location = ["UniProt loc high conf", "GO CC high conf"]

        rows = []
        for ensembl_id in summary["ensembl_gene_id"].to_list():
            target = tractability.filter((pl.col("ensembl_gene_id") == ensembl_id) & pl.col("value"))
            row = {"ensembl_gene_id": ensembl_id}
            true_labels = {}
            for modality, labels in label_order.items():
                found = target.filter(pl.col("modality") == modality)["label"].to_list()
                kept = []
                for label in labels:
                    if label in found:
                        kept.append(label)
                true_labels[modality] = kept
                row[label_columns[modality]] = "; ".join(kept)

            # small molecule: drug or ligand, else pocket only; antibody: drug, else location only; oligo: drug only
            modalities = []
            sm_strong = False
            sm_predicted = False
            for label in true_labels["SM"]:
                if label in sm_tractable:
                    sm_strong = True
                else:
                    sm_predicted = True
            if sm_strong:
                modalities.append("Small molecule")
            elif sm_predicted:
                modalities.append("Small molecule (pocket only)")
            ab_clinical = False
            ab_located = False
            for label in true_labels["AB"]:
                if label in clinical:
                    ab_clinical = True
                if label in ab_location:
                    ab_located = True
            if ab_clinical:
                modalities.append("Antibody")
            elif ab_located:
                modalities.append("Antibody (location only)")
            if len(true_labels["OC"]) > 0:
                modalities.append("Oligonucleotide")
            row["tractable_modalities"] = "; ".join(modalities)
            rows.append(row)

        labels_table = pl.DataFrame(rows, schema={
            "ensembl_gene_id": pl.String,
            "ot_small_molecule_labels": pl.String,
            "ot_antibody_labels": pl.String,
            "ot_oligonucleotide_labels": pl.String,
            "tractable_modalities": pl.String,
        })
        summary = summary.join(labels_table, on="ensembl_gene_id", how="left")
        return summary


class TargetDruggability:

    """
    Runs every source for a list of Ensembl IDs and returns dataframes (no files written):
    - target_druggability: one row per input target
    - target_drugs: one row per target x drug
    - sources: release of each source
    - raw: the untouched tables from each source
    """

    def __init__(self):
        self.dg = DruggableGenome()
        self.kd = KnownDrugs()
        self.tr = Tractability()

    def run(self, ensembl_ids, pchembl_thresh=6):

        # 1. Pharos: the base table, every input target, carries Ensembl <-> UniProt
        ensembl_ids = list(dict.fromkeys(ensembl_ids)) # drop repeats, keep order
        base = self.dg.run(ensembl_ids)
        base = base.rename({"sym": "symbol", "fam": "family"})

        # 2. ChEMBL is queried by UniProt and mapped back to Ensembl through base
        uniprots = base["uniprot"].drop_nulls().unique().to_list()
        chembl = self.kd.run(uniprots, pchembl_thresh=pchembl_thresh)
        chembl_summary = self.kd.summarise_targets(chembl)
        drugs = self.kd.drug_table(chembl)

        # 3. Open Targets is queried by Ensembl directly
        tractability = self.tr.fetch(ensembl_ids)
        ot_summary = self.tr.summarise(tractability)

        # 4. target table: left join onto base so no input target is lost
        target_table = base.select(["ensembl_gene_id", "symbol", "uniprot", "tdl", "family", "pharos_multiple_matches", "not_in_pharos"])
        target_table = target_table.join(chembl_summary, left_on="uniprot", right_on="accession", how="left")
        target_table = target_table.join(ot_summary, on="ensembl_gene_id", how="left")
        count_columns = ["n_drugs", "n_drugs_approved", "n_drugs_inhibit", "n_drugs_activate", "n_potent_compounds"]
        for column in count_columns:
            target_table = target_table.with_columns(pl.col(column).fill_null(0))

        # 5. drug table: put Ensembl ID and symbol in front
        keys = base.select(["ensembl_gene_id", "symbol", "uniprot"]).drop_nulls("uniprot")
        drug_table = keys.join(drugs, left_on="uniprot", right_on="accession", how="inner")

        # 6. releases
        today = date.today().isoformat()
        sources = pl.DataFrame({
            "source": ["Pharos/TCRD", "ChEMBL", "Open Targets Platform", "B3DB"],
            "release": [self.dg.database, self.kd.release, self.tr.release, self.kd.b3db_release],
            "accessed": [today, today, today, today],
        })

        # untouched source tables, for caching and supplementary data
        raw = {
            "pharos": base,
            "ot_tractability": tractability,
        }
        for name, table in chembl.items():
            raw["chembl_" + name] = table

        results = {
            "target_druggability": target_table,
            "target_drugs": drug_table,
            "sources": sources,
            "raw": raw,
        }

        return results


# if __name__ == "__main__":
# 
#     # Ensembl IDs: EGFR, CTSH, GPRC5B
#     td = TargetDruggability()
#     results = td.run(["ENSG00000146648", "ENSG00000103811", "ENSG00000167191"])
# 
#     print(results["target_druggability"])
#     print(results["target_drugs"])
#     print(results["sources"])
#     for name, table in results["raw"].items():
#         print(name, table.shape)
