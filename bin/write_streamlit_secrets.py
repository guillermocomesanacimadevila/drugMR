from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ENV_PATH = PROJECT_ROOT / ".env"
SECRETS_PATH = PROJECT_ROOT / ".streamlit" / "secrets.toml"

def load_env(path: Path) -> dict:
    env = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, _, value = line.partition("=")
        env[key.strip()] = value.strip()
    return env 


def write_secrets(env: dict):
    secrets = f"""[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "{env['POSTGRES_PORT']}"
database = "{env['POSTGRES_DB']}"
username = "{env['POSTGRES_USER']}"
password = "{env['POSTGRES_PASSWORD']}"
"""
    SECRETS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SECRETS_PATH.write_text(secrets, encoding="utf-8")
    print(f"[DONE] Wrote {SECRETS_PATH} from {ENV_PATH}")


def main():
    write_secrets(load_env(ENV_PATH))


if __name__ == "__main__":
    main()