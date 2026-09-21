#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 || ($1 != write && $1 != check) ]]; then
    echo 'Usage: configure-famdb.sh write|check GENE_MANIFEST DFAM_DIRECTORY' >&2
    exit 2
fi

mode=$1
manifest=$(realpath "$2")
data_dir=$(realpath "$3")
prefix="$(dirname "$manifest")/.pixi/envs/repeatmasker"

[[ -s "$manifest" ]]
[[ -d "$data_dir" ]]
mapfile -t configs < <(
    find "$prefix/share" -maxdepth 2 -type f -path '*/famdb-*/famdb.conf' -print
)
if [[ ${#configs[@]} -ne 1 ]]; then
    echo "Expected one installed famdb.conf under $prefix/share; found ${#configs[@]}." >&2
    exit 2
fi
config=${configs[0]}

if [[ $mode == check ]]; then
    configured=$(
        awk -F '=' '
            /^[[:space:]]*FAMDB_DATA_DIR[[:space:]]*=/ {
                value = $2
                sub(/^[[:space:]]+/, "", value)
                sub(/[[:space:]]+$/, "", value)
                print value
            }
        ' "$config"
    )
    [[ $configured == "$data_dir" ]] || {
        echo "famdb.conf does not bind the registered Dfam directory: $config" >&2
        exit 2
    }
    exit 0
fi

temporary=$(mktemp "${config}.tmp.XXXXXX")
trap 'rm -f "$temporary"' EXIT
printf '%s\n' \
    '# famdb.conf - managed by protist-meta-nf' \
    '# Rerun setup after every Pixi reinstall of the repeatmasker environment.' \
    '[famdb]' \
    "FAMDB_DATA_DIR = $data_dir" > "$temporary"
mv "$temporary" "$config"
trap - EXIT
