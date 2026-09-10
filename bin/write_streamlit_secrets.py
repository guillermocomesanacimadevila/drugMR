from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SECRETS_PATH = PROJECT_ROOT / ".streamlit" / "secrets.toml"
COMPOSE_PATH = PROJECT_ROOT / "docker-compose.yml"


def read_postgres_creds_from_compose():
    # docker-compose.yml is canonical - read straight from it instead of
    # hand-duplicating the same creds here, so a credential change only
    # ever has to be made in one place
    compose = yaml.safe_load(COMPOSE_PATH.read_text(encoding="utf-8"))
    env = compose["services"]["postgres"]["environment"]
    host_port = compose["services"]["postgres"]["ports"][0].split(":")[0]
    return {
        "user": str(env["POSTGRES_USER"]),
        "password": str(env["POSTGRES_PASSWORD"]),
        "database": str(env["POSTGRES_DB"]),
        "port": str(host_port),
    }


def write_secrets():
    creds = read_postgres_creds_from_compose()
    secrets = f"""[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "{creds['port']}"
database = "{creds['database']}"
username = "{creds['user']}"
password = "{creds['password']}"
"""
    SECRETS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SECRETS_PATH.write_text(secrets, encoding="utf-8")
    print(f"[DONE] Wrote {SECRETS_PATH}")


if __name__ == "__main__":
    write_secrets()
