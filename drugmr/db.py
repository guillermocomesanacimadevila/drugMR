import tomllib
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine


SECRETS_PATH = Path(__file__).resolve().parents[1] / ".streamlit" / "secrets.toml"

def load_credentials() -> dict:

    """
    Reads secrets.toml from .streamlit and 
    returns each key as a dict with respective val
    """
    
    with open(SECRETS_PATH, "rb") as f:
        secrets = tomllib.load(f)
    return secrets["connections"]["postgresql"]


def get_connections_url(db_id: str = None) -> dict:

    """
    Link credentials from streamlit app to posgtesql database
    by crafting a psql connection url from open source creds
    """

    creds = load_credentials()
    database = db_id or creds["database"]
    auth = f"{creds['username']}:{creds['password']}@" if creds["username"] else ""
    return f"postgresql://{auth}{creds['host']}:{creds['port']}/{database}"


def get_engine(db_id: str = None) -> Engine:

    """
    Get SQL engine from postgres url
    """

    return create_engine(get_connections_url(db_id=db_id))