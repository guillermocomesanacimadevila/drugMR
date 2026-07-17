import stringdb

# Map BLNK to its STRING ID
ids = stringdb.get_string_ids(["BLNK"])
blnk_id = ids.loc[0, "stringId"]

# Get BLNK's interaction partners
partners = stringdb.get_interaction_partners([blnk_id])

# BLNK + all returned partner IDs
network_ids = [blnk_id, *partners["stringId_B"].drop_duplicates().tolist()]

# Retrieve every interaction among those proteins
network = stringdb.get_network(network_ids)

print(network)
print()
print(f"Nodes: {len(set(network['stringId_A']) | set(network['stringId_B']))}")
print(f"Edges: {len(network)}")

# BRANCH INTO TWO DIFFERNT BRANCHES OF TARGATS WHICH SURVIVED PASSED SMR AND COLOC
# TREAT SEPARATELY 
# AD ONLY HITS (Y PROTEINS)
# AD + MEDIATOR HITS (X PROTEINS)
# --------
# FOR EACH SET - SMACK INTO STRINGDB
# CONFIDENCE SCORE FOR PPI > 0.40
# FOR THOSE - LEIDEN CLUSTERING - INTO CATEGORIES (CLUSTERS OF TARGETS)
# SUBMIT TO ENRICHR
# MAKE A FIGURE WITH EACH CLUSTER - EFFECT DIRECTION - KEY PATH - ETC...
# THIS DOES NOT IMPACT WHETHER WE FILTER OUT PROTEINS OR NOT 

# --------
# WITH PROTEINS WHICH SURPASSED -> COLOC/MOLOC + cis-MR
# SMR with single-cell datasets -> (SingleBrain) - largest one to date (4 differnt cohorts)
# SMR AROSS ALL CELL-TYPES -> FOR P_SMR (FDR) + HEIDI PASS ON ANY GIVEN PROTEIN:
# RUN MOLOC WITH TRAITS INVOLVED (e.g. pQTL / eQTL / GWAS/x1/x2)
