# Running as a systemd Service

For long-running deployments (lab workstation, dedicated server) we ship
a **systemd user unit** that wraps `omnipath-utils serve` so the API
restarts on failure and starts automatically at boot. This is what we
use for the public service at `utils.omnipathdb.org` and is also the
recommended setup for self-hosted instances with custom database
builds.

The unit file is part of the repository at
[`deploy/systemd/omnipath-utils.service`][unit]. It assumes:

- the package is installed into a virtualenv at `~/deploy/omnipath-utils/.venv/`,
- `~/deploy/omnipath-utils/.env` holds the runtime configuration (DB URL,
  host, port),
- a PostgreSQL container can be started via `docker compose up -d db`
  from `~/deploy/omnipath-utils/`.

[unit]: https://github.com/saezlab/omnipath-utils/blob/main/deploy/systemd/omnipath-utils.service

## One-time setup

Lay out the deploy directory the unit expects, then build the database
(see [Database Build](database.md) for build options including
multi-organism and custom resource sets):

```bash
mkdir -p ~/deploy/omnipath-utils && cd ~/deploy/omnipath-utils

# Postgres for the API to talk to
cat > docker-compose.yml << 'YAML'
services:
  db:
    image: postgres:16
    environment:
      POSTGRES_DB: omnipath_utils
      POSTGRES_USER: omnipath
      POSTGRES_PASSWORD: CHANGE_ME
    ports: ["127.0.0.1:5434:5432"]
    volumes: [pgdata:/var/lib/postgresql/data]
    restart: unless-stopped
volumes: { pgdata: {} }
YAML
docker compose up -d db

# Runtime config (read by the systemd unit)
cat > .env << 'ENV'
OMNIPATH_UTILS_DB_URL=postgresql+psycopg://omnipath:CHANGE_ME@localhost:5434/omnipath_utils
OMNIPATH_UTILS_HOST=127.0.0.1
OMNIPATH_UTILS_PORT=8083
ENV

# Virtualenv + database build
uv venv .venv --python 3.13
.venv/bin/pip install "omnipath-utils[server]"
.venv/bin/omnipath-utils build --db-url "$(grep DB_URL .env | cut -d= -f2-)" --organisms 9606
```

## Install the unit

Run as the user that will own the service (does **not** need to be
root):

```bash
git clone https://github.com/saezlab/omnipath-utils.git ~/dev/omnipath-utils
mkdir -p ~/.config/systemd/user
ln -sf ~/dev/omnipath-utils/deploy/systemd/omnipath-utils.service \
       ~/.config/systemd/user/omnipath-utils.service
systemctl --user daemon-reload
systemctl --user enable --now omnipath-utils.service
```

To survive reboots without an active login session, enable lingering
once as root:

```bash
sudo loginctl enable-linger "$USER"
```

## Custom databases

`OMNIPATH_UTILS_DB_URL` can point at any PostgreSQL instance you
control &mdash; including a database you built with a custom organism
set, custom resources, or schema overrides. The service is stateless
beyond that connection, so swapping in a new DB is a `systemctl --user
restart omnipath-utils` away.

If your environment needs additional setup (e.g. activating a Conda
profile or exporting an `LD_LIBRARY_PATH` on NixOS), drop the export
statements into `~/dev/.envrc`. The unit sources that file when
present and ignores it otherwise.

## Operate

```bash
systemctl --user status   omnipath-utils
systemctl --user restart  omnipath-utils
systemctl --user stop     omnipath-utils
journalctl  --user -u omnipath-utils -f       # live journal
tail -f ~/deploy/omnipath-utils/server.log    # uvicorn access log
```

The service binds to `127.0.0.1` by default. To expose it publicly,
front it with nginx, Caddy, or an SSH reverse tunnel and terminate TLS
there &mdash; do not bind `0.0.0.0` directly on an untrusted network.

## System-wide alternative

If you prefer a system service (e.g. a dedicated `omnipath` system
account with no login shell), copy the same unit to
`/etc/systemd/system/omnipath-utils.service`, replace `%h` with the
target home directory, add `User=` / `Group=` lines under `[Service]`,
and run `systemctl daemon-reload && systemctl enable --now
omnipath-utils.service` as root. No other changes are required.

## Multiple instances on one host

For shared development hosts where several instances of the service
need to run side-by-side (e.g. a `staging` deployment alongside a
production deployment, or a per-feature dev URL), a templated unit is
shipped at [`deploy/systemd/omnipath-utils@.service`][template].

Lay each instance out at `~/instances/<name>/` &mdash; with `src/`,
`.venv/`, `.env`, and `docker-compose.yml` (Postgres on a dedicated
port) &mdash; then enable:

```bash
systemctl --user enable --now omnipath-utils@staging.service
```

Each instance gets its own database, port, and `server.log`. The same
`loginctl enable-linger` once-per-user setup applies.

[template]: https://github.com/saezlab/omnipath-utils/blob/main/deploy/systemd/omnipath-utils@.service
