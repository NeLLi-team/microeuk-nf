#!/usr/bin/env bash
set -euo pipefail

prefix="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
share_dir="${prefix}/share/braker3-longrna"
export PATH="${share_dir}/ETP/tools:${PATH}"

exec perl "${share_dir}/ETP/bin/gmetp.pl" "$@"

