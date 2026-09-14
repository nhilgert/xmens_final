#!/usr/bin/env bash
set -euo pipefail

readonly REPO=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final
readonly BACKUP=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/_promotion_backup_20260811T231844Z
readonly MANIFEST="$BACKUP/BACKUP_MANIFEST.json"
readonly MANIFEST_SHA256=5c900eefa81c5f88ad593d0adf7e51388d99097e04ff6476c06a72f8198783c2
readonly PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python

test -t 0
test -t 1
test "$(id -u)" = 60528
cd "$REPO"
test "$(pwd -P)" = "$REPO"
[[ "$MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]]
test "$(<"$BACKUP/PREPARED_NOT_PROMOTED")" = PREPARED_NOT_PROMOTED
test ! -e "$BACKUP/PREPARATION_INCOMPLETE"
test ! -e "$BACKUP/PROMOTION_COMMITTED"
test ! -e "$BACKUP/PROMOTION_RESULT.json"

printf '%s\n' \
  '22a7ae02f71b1f8580f516763116fe668f35492f17beb596b3c31a273eb71c51  tools/apply_prepared_promotion.py' \
  '0426845841338c7c9a9f64b8b5c2fa8aac0e90eca7314d3d86968786f4599777  tools/recover_prepared_promotion.py' |
  sha256sum --check --strict -
printf '%s  %s\n' "$MANIFEST_SHA256" "$MANIFEST" |
  sha256sum --check --strict -

echo 'Validated prepared v2 receipt. No downstream job will be launched.'
exec srun --job-name=xmens_promote_v2 --partition=eddy \
  --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506 \
  --cpus-per-task=1 --mem=8G --time=04:00:00 --chdir="$REPO" --pty \
  /bin/bash -c '
set -euo pipefail
readonly REPO=$1
readonly BACKUP=$2
readonly MANIFEST_SHA256=$3
readonly PY=$4
readonly MANIFEST="$BACKUP/BACKUP_MANIFEST.json"
test -t 0
test -t 1
cd "$REPO"
test "$(pwd -P)" = "$REPO"
test "$(<"$BACKUP/PREPARED_NOT_PROMOTED")" = PREPARED_NOT_PROMOTED
test ! -e "$BACKUP/PREPARATION_INCOMPLETE"
test ! -e "$BACKUP/PROMOTION_COMMITTED"
test ! -e "$BACKUP/PROMOTION_RESULT.json"
printf "%s\n" \
  "22a7ae02f71b1f8580f516763116fe668f35492f17beb596b3c31a273eb71c51  tools/apply_prepared_promotion.py" \
  "0426845841338c7c9a9f64b8b5c2fa8aac0e90eca7314d3d86968786f4599777  tools/recover_prepared_promotion.py" |
  sha256sum --check --strict -
printf "%s  %s\n" "$MANIFEST_SHA256" "$MANIFEST" |
  sha256sum --check --strict -
exec "$PY" -I -B "$REPO/tools/apply_prepared_promotion.py" \
  --repo "$REPO" --manifest "$MANIFEST" \
  --expected-manifest-sha256 "$MANIFEST_SHA256"
' _ "$REPO" "$BACKUP" "$MANIFEST_SHA256" "$PY"
