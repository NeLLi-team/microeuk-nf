#!/usr/bin/env bash
set -euo pipefail

prefix="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
share_dir="${prefix}/share/braker3-longrna"
export GENEMARK_PATH="${share_dir}/ETP/bin"
export PROTHINT_PATH="${share_dir}/ETP/bin/gmes/ProtHint/bin"
export PATH="${share_dir}/ETP/tools:${PATH}"

exec perl "${share_dir}/BRAKER/scripts/braker.pl" "$@"

