#!/usr/bin/env bash
# Promote a build's Postgres data VOLUME to a present (web) instance WITHOUT
# pg_dump/restore.
#
# Build and web Postgres are the SAME image (omnipath-build-postgres:18 = the
# present image; both carry pg_roaringbitmap + rdkit). The only differences are
# runtime settings (compose `command:` flags), port, network — and which data
# volume is mounted. The data lives in a docker named volume, and those settings
# are startup flags, not baked into the data dir. So promotion = hand the volume
# to the web stack (which mounts it with web-tuned settings) — seconds, no copy.
#
#   promote-build.sh [--copy] <source-data-volume> <target-instance>
#
# default (hand-off): re-point the target stack's POSTGRES_DATA_VOLUME_NAME to
#   <source-data-volume> + restart it. Instant; CONSUMES the source build slot.
# --copy: raw-copy the source volume to a dated `promoted-*` volume first
#   (preserves the source build slot), then hand off the copy. ~minutes for a
#   raw file copy — still far faster than dump/restore (no index rebuild).
#
# One Postgres per data dir: any container using the source volume is stopped
# for the operation (and restarted afterwards when --copy).
set -euo pipefail
COPY=0; [ "${1:-}" = "--copy" ] && { COPY=1; shift; }
SRC="${1:?usage: promote-build.sh [--copy] <source-data-volume> <target-instance>}"
TGT="${2:?target present instance, e.g. dev3}"
ENVF="$HOME/instances/$TGT/.env"
export XDG_RUNTIME_DIR="/run/user/$(id -u)"
[ -f "$ENVF" ] || { echo "missing env file: $ENVF" >&2; exit 1; }
docker volume inspect "$SRC" >/dev/null 2>&1 || { echo "no such volume: $SRC" >&2; exit 1; }

# release the source volume (never two postgres on one data dir)
STOPPED=$(docker ps -q --filter volume="$SRC" || true)
for c in $STOPPED; do echo "stopping $(docker inspect -f '{{.Name}}' "$c") (uses $SRC)"; docker stop "$c" >/dev/null; done

VOL="$SRC"
if [ "$COPY" = 1 ]; then
  VOL="promoted-${TGT}-$(date +%Y%m%d-%H%M%S)"
  echo "raw-copying $SRC -> $VOL ..."
  docker volume create "$VOL" >/dev/null
  docker run --rm -v "$SRC":/from:ro -v "$VOL":/to alpine sh -c 'cp -a /from/. /to/'
  for c in $STOPPED; do docker start "$c" >/dev/null; echo "resumed source $(docker inspect -f '{{.Name}}' "$c")"; done
fi

sed -i "s|^POSTGRES_DATA_VOLUME_NAME=.*|POSTGRES_DATA_VOLUME_NAME=${VOL}|" "$ENVF"
echo "set POSTGRES_DATA_VOLUME_NAME=${VOL} in ${ENVF}"
systemctl --user restart "omnipath-present@${TGT}"
echo "restarted omnipath-present@${TGT}; verify http://127.0.0.1:\$API_PORT/resources"
