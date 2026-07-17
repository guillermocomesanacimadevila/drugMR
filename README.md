drugMR: A Multi-Fluid Multi-Omics Drug Target Discovery Pipeline

End-to-end local + cloud hybrid pipeline leveraging plasma and CSF protein QTLs (pQTLs; >10,000 proteins) from UKBB-PPP, deCODE and WU (Olink + SomaScan) for any phenotype of interest (demonstrated in Alzheimer's disease), with or without mediating biomarkers. From raw GWAS and pQTL data to single-cell target prioritisation, mediation analysis and drug safety assessment. Outputs are compiled into a production-ready dashboard.
---

## Clone the repo!

```bash
git clone https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR/
```

## Synapse configuration

Create the Synapse config file:

```bash
nano ~/.synapseConfig
```

Populate it with your Synapse info!

```bash
[default]
username = your_email@example.com
authtoken = YOUR_PERSONAL_ACCESS_TOKEN

[cache]
location = ~/.synapseCache
```

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
## Configure passwordless SSH access to Falcon

Generate an SSH key (if you do not already have one):

```bash
ssh-keygen -t ed25519 -C "drugMR"
```

Display your public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the public key to Falcon:

```bash
ssh-copy-id c.<username>@falconlogin.cf.ac.uk
```

Test the connection:

```bash
ssh c.<username>@falconlogin.cf.ac.uk
```


## Docker 

**This is simply a dev note, do NOT worry about it**

Pull the latest `drugMR` image from GHCR:

```bash
docker pull ghcr.io/guillermocomesanacimadevila/drugmr:latest
```

## Authors

**Guillermo Comesaña Cimadevila**<sup>1,2</sup>, **Marie-Joe Dib**<sup>3</sup>, **Dervis Salih**<sup>4</sup>, **Nicholas J. Bray**<sup>2</sup>, **Emily Simmonds**<sup>1</sup>, **Valentina Escott-Price**<sup>1,2</sup>

- <sup>1</sup> UK Dementia Research Institute at Cardiff University, Cardiff, UK.
- <sup>2</sup> MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University, Cardiff, UK.
- <sup>3</sup> Nascent Studio Ltd, London, UK.
- <sup>4</sup> UK Dementia Research Institute at University College London, London, UK.
