# cis-MR pipeline

Plasma proteins and CSF pQTLs from UKBB-PPP, NIAGADS + other sources for Alzheimer's disease drug discovery.

End-to-end pipeline from raw GWAS and pQTL data to a production-ready dashboard.

## Streamlit configuration

Create the Streamlit secrets file:

```bash
nano .streamlit/secrets.toml
```

Populate `.streamlit/secrets.toml` as follows:

```toml
[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "5432"
database = "xxx"
username = "xxx"
password = "xxx"
```

## Docker

Pull the latest `drugMR` image from GHCR:

```bash
docker pull ghcr.io/guillermocomesanacimadevila/drugmr:latest
```

## Falcon HPC authentication

Generate an SSH key:

```bash
ssh-keygen -t ed25519 -C "drugMR-falcon"
```

Display the public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the public key to Falcon:

```bash
ssh-copy-id c.<username>@falconlogin.cf.ac.uk
```

## Authors

Guillermo Comesaña Cimadevila, Marie-Joe Dib, Dervis Salih, Nicholas J. Bray, Emily Simmonds, Valentina Escott-Price
