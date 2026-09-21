#!/usr/bin/env bash
set -euo pipefail

: "${PIXI_PROJECT_ROOT:?Run this check through pixi run.}"
: "${CONDA_PREFIX:?Pixi did not activate an environment.}"
for tool in "$@"; do
    if ! executable=$(command -v "$tool"); then
        printf 'Missing Pixi dependency: %s\n' "$tool" >&2
        exit 1
    fi
    case "$executable" in
        "$CONDA_PREFIX"/*) printf '%s\t%s\n' "$tool" "$executable" ;;
        *) printf 'Dependency %s resolved outside the Pixi environment: %s\n' "$tool" "$executable" >&2; exit 1 ;;
    esac
done
