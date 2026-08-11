#!/usr/bin/env bash
# Review-only entry point for the manifest-bound xmens.annotation.v1 promotion.
# Default invocation is a read-only preflight. No downstream job is ever run.

set -euo pipefail
IFS=$'\n\t'
umask 0022

REPO="/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"
PYTHON="${REPO}/envs/samap/.venv/bin/python"

if [[ "$(pwd -P)" != "${REPO}" ]]; then
    printf 'ERROR: run from exactly %s (current: %s)\n' "${REPO}" "$(pwd -P)" >&2
    exit 2
fi

if [[ $# -eq 0 ]]; then
    set -- --dry-run
fi

exec "${PYTHON}" "${REPO}/tools/migration_promotion.py" promote "$@"
