#!/usr/bin/env bash
# Simple, human-run annotation promotion. This does not launch Snakemake or downstream jobs.
set -Eeuo pipefail
IFS=$'\n\t'
umask 0022
REPO="/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"
PLAN_REL="tools/promotion_plan_xmens_annotation_v1_20260810.json"
PLAN="$REPO/$PLAN_REL"
SCRIPT_REL="tools/promote_simple.sh"
MCK="/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data"
BASE_HEAD="2cdc6a50ef2c84c9c3340a86fafff94e95929fa6"
PLAN_FILE_SHA="bcfc92d45cd46e944a608ce9c20bedecf020cf39818d507506fb1d5cf92edc27"
PLAN_SHA="5c0dda6041034be2a6d941efe2430b028d143eaa0b1e40b055164972d0503ad3"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
BACKUP="$REPO/data/interim/_promotion_backup_$STAMP"
STARTED=0
die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
sha() { sha256sum -- "$1" | awk '{print $1}'; }
trap 'rc=$?; ((rc == 0)) || printf "STOPPED (moves_started=%s); backups: %s\n" "$STARTED" "$BACKUP" >&2' EXIT
[[ "$(pwd -P)" == "$REPO" ]] || die "run from exactly $REPO"
[[ "$(readlink -f -- "$0")" == "$REPO/$SCRIPT_REL" ]] || die "run the repository script directly"
[[ "${1:-}" == "PROMOTE:xmens_annotation_v1_20260810" && $# -eq 1 ]] ||
    die "usage: $SCRIPT_REL PROMOTE:xmens_annotation_v1_20260810"
for command in git jq sha256sum cp mv flock; do command -v "$command" >/dev/null || die "missing $command"; done
[[ "$(git rev-parse --show-toplevel)" == "$REPO" ]] || die "wrong Git repository"
git merge-base --is-ancestor "$BASE_HEAD" HEAD || die "reviewed base is not an ancestor of HEAD"
[[ -z "$(git status --porcelain --untracked-files=all -- "$PLAN_REL" "$SCRIPT_REL")" ]] ||
    die "plan and script must be committed and clean"
[[ "$(sha "$PLAN")" == "$PLAN_FILE_SHA" ]] || die "promotion plan file changed"
jq -e --arg p "$PLAN_SHA" --arg h "$BASE_HEAD" '.canonical_plan_sha256==$p and
  .reviewed_repository_head==$h and (.plan.mappings|length)==21 and
  ([.plan.mappings[].order]==[range(1;22)])' "$PLAN" >/dev/null || die "promotion plan contract failed"
exec 9>"$REPO/data/interim/.annotation_promotion.lock"
flock -n 9 || die "another promotion process holds the lock"
declare -a SRC DST OLD NEW BAK TMP
mapfile -t rows < <(jq -r '.plan.mappings[] | [.stage,.source.path,.source.sha256,.destination.path,.destination.old.sha256] | @tsv' "$PLAN")
for row in "${rows[@]}"; do
    IFS=$'\t' read -r stage source new destination old <<<"$row"
    SRC+=("$REPO/$source"); DST+=("$REPO/$destination"); NEW+=("$new"); OLD+=("$old")
    if [[ "$stage" == deposit ]]; then
        SRC+=("$REPO/$source"); DST+=("$MCK/${destination##*/}"); NEW+=("$new"); OLD+=("$old")
    fi
done
[[ ${#DST[@]} -eq 25 ]] || die "expected 21 reviewed targets plus four McKinley mirrors"
for i in "${!DST[@]}"; do
    [[ -f "${SRC[i]}" && ! -L "${SRC[i]}" ]] || die "bad staged source: ${SRC[i]}"
    [[ -f "${DST[i]}" && ! -L "${DST[i]}" ]] || die "bad live destination: ${DST[i]}"
    [[ "$(stat -c %h -- "${DST[i]}")" == 1 ]] || die "hard-linked destination: ${DST[i]}"
    [[ "$(stat -c %u -- "${DST[i]}")" == "$(id -u)" ]] || die "destination not owned by caller"
    [[ "$(sha "${SRC[i]}")" == "${NEW[i]}" ]] || die "staged hash changed: ${SRC[i]}"
    [[ "$(sha "${DST[i]}")" == "${OLD[i]}" ]] || die "live hash changed: ${DST[i]}"
done
mkdir -p -- "$(dirname "$BACKUP")"; mkdir -- "$BACKUP"
printf 'index\tsource\tdestination\tbackup\told_sha256\tnew_sha256\n' >"$BACKUP/MANIFEST.tsv"
for i in "${!DST[@]}"; do
    printf -v tag '%02d_%s' "$((i+1))" "$(basename "${DST[i]}")"
    BAK[i]="$BACKUP/$tag.old"; TMP[i]="${DST[i]}.promote-simple.$STAMP.tmp"
    [[ ! -e "${TMP[i]}" ]] || die "temporary path exists: ${TMP[i]}"
    cp --reflink=never --sparse=always --preserve=all -- "${DST[i]}" "${BAK[i]}"
    cp --reflink=never --sparse=always --preserve=mode,timestamps -- "${SRC[i]}" "${TMP[i]}"
    chmod --reference="${DST[i]}" -- "${TMP[i]}"; chgrp --reference="${DST[i]}" -- "${TMP[i]}"
    [[ "$(sha "${BAK[i]}")" == "${OLD[i]}" ]] || die "backup verification failed"
    [[ "$(sha "${TMP[i]}")" == "${NEW[i]}" ]] || die "temp verification failed"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$((i+1))" "${SRC[i]}" "${DST[i]}" \
        "${BAK[i]}" "${OLD[i]}" "${NEW[i]}" >>"$BACKUP/MANIFEST.tsv"
done
for i in "${!DST[@]}"; do
    [[ "$(sha "${BAK[i]}")" == "${OLD[i]}" && "$(sha "${TMP[i]}")" == "${NEW[i]}" ]] ||
        die "prepared copy changed: ${DST[i]}"
done
printf 'Prepared and verified 25 backups and 25 adjacent install files; beginning atomic moves.\n'
STARTED=1
for i in "${!DST[@]}"; do
    [[ "$(sha "${DST[i]}")" == "${OLD[i]}" ]] || die "destination changed before move: ${DST[i]}"
    mv -fT -- "${TMP[i]}" "${DST[i]}"
    [[ "$(sha "${DST[i]}")" == "${NEW[i]}" ]] || die "installed hash failed: ${DST[i]}"
done
for i in "${!DST[@]}"; do [[ "$(sha "${DST[i]}")" == "${NEW[i]}" ]] || die "final hash failed"; done
printf 'PASS: 21 reviewed live targets plus four independent McKinley mirrors promoted.\n'
printf 'Verified backups and manifest: %s\nNo downstream job was run.\n' "$BACKUP"
