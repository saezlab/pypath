#!/usr/bin/env bash
# Start the local Parquet API and explorer. No resource or database builds run.
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_dir"
api_port="${1:-8085}"
web_port="${2:-5173}"
data_root="${OMNIPATH_DATA_ROOT:-data}"

uv run --frozen --package omnipath-api omnipath-api --host 127.0.0.1 --port "$api_port" --data-root "$data_root" &
api_pid=$!
cleanup() {
    kill "$api_pid" 2>/dev/null || true
    wait "$api_pid" 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

API_SERVICE_URL="http://127.0.0.1:$api_port" \
    pnpm --dir packages/web dev --host 127.0.0.1 --port "$web_port"
