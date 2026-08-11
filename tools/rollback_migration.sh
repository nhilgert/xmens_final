#!/usr/bin/env bash
# Review-only entry point for restoring the 21 reviewed pre-migration objects.
# A dry-run still needs the promotion id and printed backup-manifest SHA-256.

set -euo pipefail
IFS=$'\n\t'
umask 0022

REPO="/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"
PYTHON="${REPO}/envs/samap/.venv/bin/python"

if [[ "$(pwd -P)" != "${REPO}" ]]; then
    printf 'ERROR: run from exactly %s (current: %s)\n' "${REPO}" "$(pwd -P)" >&2
    exit 2
fi

exec "${PYTHON}" "${REPO}/tools/migration_promotion.py" rollback "$@"
