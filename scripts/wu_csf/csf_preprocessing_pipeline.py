#!/urs/bin/env python3

# 1. Load WU CSF manifest
# 2. Keep:
#       - Human targets only
#       - Autosomal targets only: chromosomes 1–22
#       - Successfully mapped genes
#       - Valid chromosome, Begin and End coordinates
#       - Valid GCST accession and harmonised FTP URL
#
# 3. Create a stable assay identifier:
#       gene_uniprot_seqid
#
#    Do not use only geneID_uniprotID because multiple aptamers can measure
#    the same protein. SeqId should remain part of the identifier.
#
# 4. Define cis region:
#       cis_start = max(Begin - 1_000_000, 1)
#       cis_end   = End + 1_000_000
#
# 5. For each SeqId / GCST accession:
#       - Download the harmonised GWAS .tsv.gz
#       - Extract variants within the encoding gene cis region
#       - Standardise columns
#       - QC
#       - Save the cleaned cis region as Parquet
#       - Delete the downloaded .tsv.gz
