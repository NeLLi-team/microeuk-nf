#!/usr/bin/env bash
set -euo pipefail

prefix="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
prothint_dir="${prefix}/share/braker3-longrna/ETP/bin/gmes/ProtHint"

exec python "${prothint_dir}/bin/prothint.py" "$@"

