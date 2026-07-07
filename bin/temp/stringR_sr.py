import stringdb

# Map gene symbol -> STRING ID
ids = stringdb.get_string_ids(["ADAM10"])

string_id = ids.loc[0, "stringId"]

# Retrieve interaction partners
partners = stringdb.get_interaction_partners([string_id])

print(partners)

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
