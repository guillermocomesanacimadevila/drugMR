# cis-MR pipeline

Plasma proteins and CSF (pQTLs) from UKBB-PPP, NIAGADS + other sources for Alzheimer's disease drug discovery

- End-to-end pipeline from raw GWAS and pQTL data to a production-ready dashboard

´´´´bash
nano .streamlit/secrets.toml
´´´


# .streamlit/secrets.toml

[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "5432"
database = "xxx"
username = "xxx"
password = "xxx"y

---

## Authors

Guillermo Comesaña Cimadevila, Marie-Joe Dib, Dervis Salih, Nicholas J. Bray, Emily Simmonds, Valentina Escott-Price



## Docker

docker build -t drugmr:latest -f env/Dockerfile .
