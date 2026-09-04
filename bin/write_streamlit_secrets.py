from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SECRETS_PATH = PROJECT_ROOT / ".streamlit" / "secrets.toml"

# must match docker-compose.yml's postgres service exactly
POSTGRES_USER = "drugmr_user"
POSTGRES_PASSWORD = "drugmr_ukdri"
POSTGRES_DB = "drugmr"
POSTGRES_PORT = "5433"


def write_secrets():
    secrets = f"""[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "{POSTGRES_PORT}"
database = "{POSTGRES_DB}"
username = "{POSTGRES_USER}"
password = "{POSTGRES_PASSWORD}"
"""
    SECRETS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SECRETS_PATH.write_text(secrets, encoding="utf-8")
    print(f"[DONE] Wrote {SECRETS_PATH}")


if __name__ == "__main__":
    write_secrets()
