# AFK annotation migration regeneration and promotion preparation

Date: 2026-08-11 (America/New_York)

## Executive status

- Promotion performed: **NO**.
- SAMap run: **NO**.
- Live protected objects modified: **NO**.
- Selected promotion source: `data/interim/migrated/xmens_annotation_v2_20260811/`.
- v1 fallback: **REJECTED** because its two staged whole-atlas finals were overwritten after v1 certification while its deposits/fib/dec artifacts remained older.
- Final preparation status: **PREPARED_NOT_PROMOTED — PASS** (25/25 backups verified; human approval not run).

## Phase 1 — isolated v2 regeneration

The standalone workflow `workflow/migration/Snakefile_v2_20260811` wrote only to:

- `data/interim/migrated/xmens_annotation_v2_20260811/`
- `dataverse/objects_migrated/xmens_annotation_v2_20260811/`
- `logs/migration/xmens_annotation_v2_20260811/`

The initial dry-run scheduled exactly 16 jobs: 1 protected baseline capture, 11 lineage reannotations (04c), 2 merges (04d), and 2 whole-atlas re-embeds (04e). It scheduled zero 04b lineage refits. The 04c/04d jobs ran CPU-only; only the two 04e jobs requested one generic GPU each. Slurm allocated an A40 for mouse and an A100-SXM4-80GB for human.

### Count gates

| Species | Gate | Expected | Actual | Result |
|---|---:|---:|---:|---|
| Human | Total cells | 66,830 | 66,830 | PASS |
| Mouse | Total cells | 57,033 | 57,033 | PASS |
| Human | Decidual cells | 21,648 | 21,648 | PASS |
| Human | Endometrial fibroblasts | 25,885 | 25,885 | PASS |
| Human | Epithelial cells | 10,699 | 10,699 | PASS |
| Human | Immune cells | 4,709 | 4,709 | PASS |
| Human | Smooth-muscle cells | 2,084 | 2,084 | PASS |
| Human | Endothelial cells | 1,334 | 1,334 | PASS |
| Human | Perivascular cells | 471 | 471 | PASS |
| Mouse | Decidual cells | 22,919 | 22,919 | PASS |
| Mouse | Endometrial fibroblasts | 20,278 | 20,278 | PASS |
| Mouse | Epithelial cells | 6,265 | 6,265 | PASS |
| Mouse | Immune cells | 2,091 | 2,091 | PASS |
| Mouse | Smooth-muscle cells | 1,443 | 1,443 | PASS |
| Mouse | Endothelial cells | 3,487 | 3,487 | PASS |
| Mouse | Mesothelial cells | 550 | 550 | PASS |
| Human | Empty `cell_subtype` | 0 | 0 | PASS |
| Mouse | Empty `cell_subtype` | 0 | 0 | PASS |
| Human | Decidual coarse fill | 21,648 | 21,648 | PASS |
| Human | Fibroblast coarse fill | 25,885 | 25,885 | PASS |
| Human | Non-ciliated epithelial coarse fill | 10,608 | 10,608 | PASS |
| Human | Perivascular coarse fill | 471 | 471 | PASS |
| Human | Smooth-muscle coarse fill | 2,084 | 2,084 | PASS |
| Mouse | Epithelial coarse fill | 6,265 | 6,265 | PASS |
| Mouse | Mesothelial coarse fill | 550 | 550 | PASS |
| Mouse | Smooth-muscle coarse fill | 1,443 | 1,443 | PASS |

Full validation result: **PASS**. Both bundle validators and both annotation-contract manifests passed; all 21 unique staged objects matched their recorded SHA-256 and size pins. Human fib/dec contained 47,533 cells and mouse fib/dec contained 43,197 cells. Every PASS count gate had delta zero; the only INFO rows were the expected legacy-reference agreement comparisons.

Protected-live post-run result: **PASS** (`26/26`, `changes=[]`). An independent rehash also matched content, size, inode, and mtime for all 26 protected files and matched all 18 frozen small inputs.

After all 25 backups were prepared, a second full protected-live re-hash (job
`38455741`) completed `0:0` in 21:04 and again returned **PASS** (`26/26`,
`changes=[]`). Its report SHA-256 is
`0b6b24e3bb71ec53a109799736432c21a4c3e7e7b9ad9c8e57c7a476462d4b78`;
stderr was empty. This final gate proves the backup preparation itself did not
modify a protected live object.

Fresh whole-atlas SHA-256 values:

- Human: `1bc8b2a8bf4fafab7ee189ad9bc04157bb0e270a305e839effaa605f7cd4c1c8`
- Mouse: `c297c931320c37927550819bcf6f538b4731bdd9f17f0e81236a74d306742ca0`

Both refined-to-final comparisons preserved cell IDs, gene IDs, and canonical annotation arrays in exact order. Both 30-dimensional scVI matrices and two-dimensional UMAP matrices were finite; all 21 Leiden columns were complete.

## Phase 2 — public DRAFT notebook execution

The original 11 files in `notebooks/public/` were not edited. Byte-identical copies were executed under `data/interim/_draft_exec_inline_20260811/`, with HOME, Jupyter, matplotlib, Numba, Torch, CUDA/Triton, Python bytecode, XDG, and temporary paths isolated below that scratch root.

The embedded `ANNOTATION` dictionaries matched their committed decision CSVs in all 11 notebooks. Each draft embeds the older vocabulary digest `df2709…`, whereas the coarse-fill vocabulary is `d02b6812…`; therefore the expected honest outcome is one final-cell stale-vocabulary `AssertionError` per notebook after all six code cells and plots execute.

| Species | Lineage | Notebook status | Diagnostic classification | CSV map | Plot outputs | Key error |
|---|---|---|---|---|---:|---|
| Human | epithelial | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Human | immune | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Human | perivascular | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Human | smc | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Human | stromal | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Mouse | epithelial | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Mouse | immune | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Mouse | perimetrium | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Mouse | smc | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Mouse | stromal | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |
| Mouse | vascular | ERROR | EXPECTED_ERROR | PASS | 3 | `vocabulary changed; regenerate this draft` |

Aggregate audit: **PASS** (`integrity_status=PASS`, `diagnostic_contract=PASS`, 11 ERROR, 11 EXPECTED_ERROR, 0 clean PASS). All six code cells executed in every notebook, source cells remained unchanged, and lineage IDs/types/subtypes matched v2. The first scratch attempt used `MPLBACKEND=Agg`; it honestly failed the diagnostic contract because that backend suppressed all inline PNG payloads. A separate reviewed scratch run switched only to the headless inline backend and produced three images per notebook. Both attempts remain available for audit.

## Phase 3 — promotion preparation only

Chosen source: `data/interim/migrated/xmens_annotation_v2_20260811/`, because v2 passed the count, schema, object, export, and protected-live gates. No automatic fallback was permitted.

Backup directory: `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/_promotion_backup_20260811T231844Z`

Backup manifest SHA-256: `5c900eefa81c5f88ad593d0adf7e51388d99097e04ff6476c06a72f8198783c2`

Backup job `38428718` completed `0:0` in 1:21:02. It copied exactly
32,189,581,132 bytes of reviewed-old data. `BACKUP_MANIFEST.tsv` has SHA-256
`3b1bb96efd0a5bda713c802cf7d6388366d7d7459301e06c200d957a4b3aa0a3`;
the exact `PREPARED_NOT_PROMOTED\n` receipt has SHA-256
`3b495a93c57c29675b1f337ef475c3596fc5d3d1c02d9f8b7cdaef94bc72f33b`.
No `PREPARATION_INCOMPLETE` or promotion/recovery marker exists.
An independent read-only audit reconciled JSON to TSV, all 21 source
certifications, all 25 current destinations, and all 25 backups, including
content and recorded metadata; it found no mismatch or SUSPECT result. GPFS
reported device 63 from the compute node and device 55 from the login node,
but source/destination/backup were same-device within each host, so this is a
benign host-view difference rather than a cross-filesystem copy.

All 25 independent backups (21 reviewed live destinations plus 4 McKinley mirrors) were copied and SHA-256 verified. No live pathname was replaced. The dry-run verified exact source/destination mapping, old and new hashes, ownership/mode, same-device atomic-install compatibility, writable destination parents, and dependency order.

| Physical rows | Kind | Verification |
|---:|---|---|
| 1–11 | Reannotated lineage objects | PASS — old SHA-256 and size verified |
| 12–13 | Refined atlases | PASS — old SHA-256 and size verified |
| 14–15 | Final atlases | PASS — old SHA-256 and size verified |
| 16–17 | Fib/dec subsets | PASS — old SHA-256 and size verified |
| 18, 20, 22, 24 | Dataverse deposits | PASS — old SHA-256 and size verified |
| 19, 21, 23, 25 | Independent McKinley mirrors | PASS — old SHA-256 and size verified |

### Complete promotion dry-run

| Order | Destination | Staged source | Old SHA-256 | New SHA-256 | Result |
|---:|---|---|---|---|---|
| 01 | `data/interim/reannotated/human_stromal.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/human_stromal.h5ad` | `05b152581e474b8ecd49fbf6e5d4b34cca8684b40165e7c7020590557f9d6c52` | `0d33b6176ce9573375bd1f681a2914a5971714d40d702233781f999a455e9c40` | would promote |
| 02 | `data/interim/reannotated/human_epithelial.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/human_epithelial.h5ad` | `61566aa010011da0570ca6d931ed3949a4067a03569b2ec4e083337672383cdf` | `c3a873ca8ecff7341adec891834e0a7c2243341bba7413de5e3eea50e6a08182` | would promote |
| 03 | `data/interim/reannotated/human_immune.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/human_immune.h5ad` | `d49d1373ed92322d3e6b9af0cbaaaeec0f19838ff1fa522015c808664d9a1b07` | `810f2a391cac2d3065c06f7690841bc7b7ce9a437de75838e294c0c279f59575` | would promote |
| 04 | `data/interim/reannotated/human_perivascular.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/human_perivascular.h5ad` | `9b6f30791fdc51061433d433d78fb4612fe7920cf7380df02c465d66859895e1` | `af73161ee5c9cfc9f2cf43c38a41d366165a58a5b21a1facbcd3690c153d6dcc` | would promote |
| 05 | `data/interim/reannotated/human_smc.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/human_smc.h5ad` | `373c7a8b8c15f92244f97bc5e0b8581d0d028d62be176d5b2dd6ce578cdff379` | `d62f00723bd2a6d233aee1d4ba2d3129179ef49eb51fa5fa6e5406e723ab0d3c` | would promote |
| 06 | `data/interim/reannotated/mouse_stromal.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/mouse_stromal.h5ad` | `7016f72fe8d03142dd5c2d59fb26e71efb1af71ef04cfd02880cfc71a8c5e96b` | `496162972869d8babf72acf04b696b18bb07d5d7e5841f3678061102ff031d88` | would promote |
| 07 | `data/interim/reannotated/mouse_epithelial.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/mouse_epithelial.h5ad` | `6a507124c09bf16d406b1d6226be1a7113fe56fb076c4bc6e5f37e4930647d43` | `fe649c3a4520862bce9d240a8923d30a7b371e1a0ddf265a58dfa0ea162c4030` | would promote |
| 08 | `data/interim/reannotated/mouse_immune.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/mouse_immune.h5ad` | `fba4cdcba00dbabc71b58ddf87edc617dc3f63da7502fd4d42b5dcdb0ae197cf` | `08778b953559738f0afc23ffbe3ec3e4921fae39f7ad23e88507ec594da8cdfa` | would promote |
| 09 | `data/interim/reannotated/mouse_vascular.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/mouse_vascular.h5ad` | `3174cad1f47f2879ff6369118776053d14d860db08fcaabcd42ffe2ea0e4dff3` | `9966b9f57f0d1424189b8803ed34859444d3ed2eba2f00ae200e0888b6ca3b79` | would promote |
| 10 | `data/interim/reannotated/mouse_smc.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/mouse_smc.h5ad` | `dfbe2ec008669d9b1a3d72b7a233ac35da1163f399717c59491b2fe0222f012f` | `f8052a23a401bfc75205ada85cb74fd60ba30caf4e4a317775cae8af7b335448` | would promote |
| 11 | `data/interim/reannotated/mouse_perimetrium.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/reannotated_migrated/mouse_perimetrium.h5ad` | `efdc119ad8cd75d009b96c2ee297a8b7ef0c336d5a2069eeef830bc3ded57036` | `fc29226f7c0bc6d5a4ff20313b13cdc4d1051ced631ffb8db2c5183714cf7d12` | would promote |
| 12 | `data/interim/refined_human.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/refined_human_migrated.h5ad` | `a31f53b717083e10bccbc58595b8538f962dcfbeb235d239a175a5bb777fcc6a` | `eb2d08c910789ab23ad2812e0bd605dc1c8a5c5fbc618f63ba06c2c441e5a898` | would promote |
| 13 | `data/interim/refined_mouse.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/refined_mouse_migrated.h5ad` | `1a0b1b4394ccbb48b3fa5ef31d5afa00c99427f1d476ea3e751f39153297ab55` | `ddb28b5e5684074ebaa99b436cd917e9a70ab751b69c9700a72fed35bfdccc0b` | would promote |
| 14 | `data/interim/final_human.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/final_human_migrated.h5ad` | `0b9582d374ecdc6a955982d7edbfe81eb69d44dc26f7de1515effdc3e6da09dd` | `1bc8b2a8bf4fafab7ee189ad9bc04157bb0e270a305e839effaa605f7cd4c1c8` | would promote |
| 15 | `data/interim/final_mouse.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/final_mouse_migrated.h5ad` | `7cfa0ebda34812fbfbc9ea4fbe23888e32997ad0bc8263634abc7f9f30ad7116` | `c297c931320c37927550819bcf6f538b4731bdd9f17f0e81236a74d306742ca0` | would promote |
| 16 | `data/interim/final_human_fib_dec.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/final_human_fib_dec_migrated.h5ad` | `511172731e5f7cefc910ad4fd021fc9e8545ed9d4e639ee2bdeae282918e9af5` | `ab138d33ad82bb4139a5a5b80b72a65f28d3ff5c5030ac09cbbe70311a93e000` | would promote |
| 17 | `data/interim/final_mouse_fib_dec.h5ad` | `data/interim/migrated/xmens_annotation_v2_20260811/final_mouse_fib_dec_migrated.h5ad` | `147a9eddba4c4982b75e03d2557d33c17d95d76fee124b376c033fd2e4e11923` | `13cd2730fd792186532545e88e0c71aefe0d1ea59d0ef14a67e43637629c33ed` | would promote |
| 18 | `dataverse/objects/2026_08_05_Cevrim_Human_All_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_Human_All_Cells.h5ad` | `0b9582d374ecdc6a955982d7edbfe81eb69d44dc26f7de1515effdc3e6da09dd` | `1bc8b2a8bf4fafab7ee189ad9bc04157bb0e270a305e839effaa605f7cd4c1c8` | would promote |
| 19 | `/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/2026_08_05_Cevrim_Human_All_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_Human_All_Cells.h5ad` | `0b9582d374ecdc6a955982d7edbfe81eb69d44dc26f7de1515effdc3e6da09dd` | `1bc8b2a8bf4fafab7ee189ad9bc04157bb0e270a305e839effaa605f7cd4c1c8` | would promote |
| 20 | `dataverse/objects/2026_08_05_Cevrim_Human_Only_Fib_Dec_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_Human_Only_Fib_Dec_Cells.h5ad` | `511172731e5f7cefc910ad4fd021fc9e8545ed9d4e639ee2bdeae282918e9af5` | `ab138d33ad82bb4139a5a5b80b72a65f28d3ff5c5030ac09cbbe70311a93e000` | would promote |
| 21 | `/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/2026_08_05_Cevrim_Human_Only_Fib_Dec_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_Human_Only_Fib_Dec_Cells.h5ad` | `511172731e5f7cefc910ad4fd021fc9e8545ed9d4e639ee2bdeae282918e9af5` | `ab138d33ad82bb4139a5a5b80b72a65f28d3ff5c5030ac09cbbe70311a93e000` | would promote |
| 22 | `dataverse/objects/2026_08_05_Cevrim_XMens_All_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_XMens_All_Cells.h5ad` | `7cfa0ebda34812fbfbc9ea4fbe23888e32997ad0bc8263634abc7f9f30ad7116` | `c297c931320c37927550819bcf6f538b4731bdd9f17f0e81236a74d306742ca0` | would promote |
| 23 | `/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/2026_08_05_Cevrim_XMens_All_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_XMens_All_Cells.h5ad` | `7cfa0ebda34812fbfbc9ea4fbe23888e32997ad0bc8263634abc7f9f30ad7116` | `c297c931320c37927550819bcf6f538b4731bdd9f17f0e81236a74d306742ca0` | would promote |
| 24 | `dataverse/objects/2026_08_05_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` | `147a9eddba4c4982b75e03d2557d33c17d95d76fee124b376c033fd2e4e11923` | `13cd2730fd792186532545e88e0c71aefe0d1ea59d0ef14a67e43637629c33ed` | would promote |
| 25 | `/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/2026_08_05_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` | `dataverse/objects_migrated/xmens_annotation_v2_20260811/2026_08_05_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` | `147a9eddba4c4982b75e03d2557d33c17d95d76fee124b376c033fd2e4e11923` | `13cd2730fd792186532545e88e0c71aefe0d1ea59d0ef14a67e43637629c33ed` | would promote |

`tools/APPROVE_AND_PROMOTE.sh` was written and reviewed but **not run**. It pins the exact prepared manifest and audited helper hashes, requires an interactive typed confirmation, creates verified same-directory promotion and rollback temporaries, journals every atomic replacement, automatically rolls back ordinary failures/signals, and launches no downstream job. A separate pinned recovery helper handles a hard stop before the durable commit marker.

### Human promotion command

Run only after inspecting this report and the backup manifest:

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final
./tools/APPROVE_AND_PROMOTE.sh
```

Keep the prepared backup directory after a successful promotion until the promoted
objects and durable commit/result manifests have been reviewed.

### Recovery command for an interrupted human promotion

Run this only if the approval helper was started, `RECOVERY_REQUIRED` exists, and
`PROMOTION_COMMITTED` does not exist. It restores and verifies all 25 reviewed old
destinations from the prepared backups and requires its own typed confirmation.

```bash
set -euo pipefail
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final
test -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/_promotion_backup_20260811T231844Z/RECOVERY_REQUIRED
test ! -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/_promotion_backup_20260811T231844Z/PROMOTION_COMMITTED
srun --job-name=xmens_recover_v2 --partition=eddy \
  --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506 \
  --cpus-per-task=1 --mem=8G --time=04:00:00 \
  --chdir=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final --pty \
  /bin/bash -c '
set -euo pipefail
readonly REPO=$1 BACKUP=$2 EXPECTED=$3
cd "$REPO"
printf "%s  %s\n" \
  0426845841338c7c9a9f64b8b5c2fa8aac0e90eca7314d3d86968786f4599777 \
  tools/recover_prepared_promotion.py | sha256sum --check --strict -
printf "%s  %s\n" "$EXPECTED" "$BACKUP/BACKUP_MANIFEST.json" |
  sha256sum --check --strict -
exec /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python -I -B \
  "$REPO/tools/recover_prepared_promotion.py" \
  --repo "$REPO" --manifest "$BACKUP/BACKUP_MANIFEST.json" \
  --expected-manifest-sha256 "$EXPECTED"
' _ /n/lab_storage/eddy_lab/users/nhilgert/xmens_final \
  /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/_promotion_backup_20260811T231844Z \
  5c900eefa81c5f88ad593d0adf7e51388d99097e04ff6476c06a72f8198783c2
```

### Safe discard command

This removes only the unpromoted v2 staging, v2 deposits, migration logs, notebook-execution scratch, and the exact prepared backup directory. Do not run it after promotion.

```bash
set -euo pipefail
readonly REPO=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final
readonly BACKUP=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/_promotion_backup_20260811T231844Z
readonly EXPECTED_MANIFEST=5c900eefa81c5f88ad593d0adf7e51388d99097e04ff6476c06a72f8198783c2
readonly PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python
cd "$REPO"
test "$(pwd -P)" = "$REPO"
test "$BACKUP" = "$REPO/data/interim/_promotion_backup_20260811T231844Z"
exec 9>"$REPO/data/interim/.annotation_promotion.lock"
flock -n 9
test "$(<"$BACKUP/PREPARED_NOT_PROMOTED")" = PREPARED_NOT_PROMOTED
test ! -e "$BACKUP/PREPARATION_INCOMPLETE"
printf '%s  %s\n' "$EXPECTED_MANIFEST" "$BACKUP/BACKUP_MANIFEST.json" |
  sha256sum --check --strict -

# Refuse to discard if promotion/recovery ever started, if any transaction temp
# remains, or if any of the 25 live destinations differs from reviewed-old state.
"$PY" -I -B - "$BACKUP/BACKUP_MANIFEST.json" "$EXPECTED_MANIFEST" <<'PY'
import hashlib, json, os, stat, sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
expected_manifest = sys.argv[2]
backup = manifest_path.parent

def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest_bytes = manifest_path.read_bytes()
if hashlib.sha256(manifest_bytes).hexdigest() != expected_manifest:
    raise SystemExit("REFUSE: prepared manifest changed")
allowed_backup_entries = {
    "repo", "mckinley", "BACKUP_MANIFEST.json", "BACKUP_MANIFEST.tsv",
    "PREPARED_NOT_PROMOTED",
}
unexpected = sorted(entry.name for entry in backup.iterdir()
                    if entry.name not in allowed_backup_entries)
if unexpected:
    raise SystemExit(f"REFUSE: transaction or unexpected backup artifacts exist: {unexpected}")
manifest = json.loads(manifest_bytes)
rows = manifest.get("rows", [])
if (manifest.get("status") != "PREPARED_NOT_PROMOTED"
        or manifest.get("run_id") != "xmens_annotation_v2_20260811"
        or len(rows) != 25
        or [row.get("physical_order") for row in rows] != list(range(1, 26))):
    raise SystemExit("REFUSE: prepared manifest contract is wrong")
for row in rows:
    path = Path(row["destination"])
    st = os.lstat(path)
    if not stat.S_ISREG(st.st_mode) or st.st_nlink != 1:
        raise SystemExit(f"REFUSE: live destination is not an independent regular file: {path}")
    old = row["old_observed"]
    observed = {
        "size": st.st_size, "mode": f"{stat.S_IMODE(st.st_mode):04o}",
        "uid": st.st_uid, "gid": st.st_gid, "mtime_ns": st.st_mtime_ns,
        "inode": st.st_ino,
    }
    for key, value in observed.items():
        if old.get(key) != value:
            raise SystemExit(f"REFUSE: old metadata drift ({key}): {path}")
    if digest(path) != old["sha256"]:
        raise SystemExit(f"REFUSE: old content drift: {path}")
    for pattern in (
        f".{path.name}.promote.*.tmp", f".{path.name}.restore.*.tmp",
        f".{path.name}.manual_restore.*.tmp",
        f".{path.name}.emergency_restore.*.tmp",
    ):
        transaction_temps = list(path.parent.glob(pattern))
        if transaction_temps:
            raise SystemExit(f"REFUSE: transaction temp exists: {transaction_temps}")
print("PASS: promotion was never attempted and all 25 live destinations remain reviewed-old")
PY

# Still under the shared exclusive lock: re-scan the backup root immediately
# before deletion so no transaction evidence can appear between verification and rm.
unexpected=$(find "$BACKUP" -mindepth 1 -maxdepth 1 \
  ! -name repo ! -name mckinley \
  ! -name BACKUP_MANIFEST.json ! -name BACKUP_MANIFEST.tsv \
  ! -name PREPARED_NOT_PROMOTED -print -quit)
test -z "$unexpected"

rm -rf -- \
  "$REPO/data/interim/migrated/xmens_annotation_v2_20260811" \
  "$REPO/dataverse/objects_migrated/xmens_annotation_v2_20260811" \
  "$REPO/logs/migration/xmens_annotation_v2_20260811" \
  "$REPO/data/interim/_draft_exec_20260811" \
  "$REPO/data/interim/_draft_exec_inline_20260811" \
  "$BACKUP"
```

## Findings and caveats

- SUSPECT items: **none**. The 11 notebook errors are expected stale-vocabulary
  diagnostics, not content/count SUSPECT results.
- v1 is not a safe fallback: its certified final-human/final-mouse staged hashes no longer match the v1 promotion plan, and the v1 root is internally inconsistent with its older deposits.
- The notebook drafts are valid decision/plot scaffolds but are intentionally reported `ERROR / EXPECTED_ERROR` because their embedded vocabulary digest is stale. They should be regenerated later; the originals were left untouched.
- Promotion is per-file atomic, not a filesystem-wide 25-file atomic operation. The approval helper therefore pre-stages complete promotion and restoration sets, uses a durable recovery marker/journal, and blocks/rolls back on any pre-commit failure.
- No SAMap or downstream analysis was run.
- Inventoried and removed three stale home Snakemake runtime-cache trees:
  `tmp6480245y` (60,188 bytes), `tmpc8g830mr` (59,501 bytes), and
  `tmpq9dw19y1` (12,693 bytes). The two August trees came from earlier
  unpinned Snakemake help/version probes; all production workflow and notebook
  jobs in this task used caches below lab storage. A post-cleanup search found
  no remaining `tmp*` tree at that runtime-cache location.
- No stale `/n/eddy_lab/` path was used by the production run; every production
  repo, data, environment, cache, log, and scratch path resolved under
  `/n/lab_storage/eddy_lab/`. The two earlier help/version-probe cache writes are
  the explicitly inventoried and removed exception described above.

## Reproducibility pins

- v2 workflow SHA-256: `3195d68cf4066247d75ff0c9bbeda21f763862854c11a9519b391c1d59ac1106`
- annotation producer SHA-256: `e57cd6e6b66ece138e770048060e7eece2eb6020395bc21b73aab4e78812b480`
- vocabulary SHA-256: `d02b68120a62d1764b6f0473976fcd72a85f5c30a081b3a857db9215e6a38847`
- draft prepare/runner/audit SHA-256: `c75734d201ede6370fa33135b962e499092a0fb5ec670a7f351dc9d12277372f`, `7ce23b2142148749e974202974d4f6f0b9a987f94d8f4bced94c4a2332139d84`, `4bc17f2b83ac5aaf8a80eeea3310a7ea63e668ae991f2577af680a769d4c49fd`
- corrected inline draft prepare/runner/audit SHA-256: `305c5439e564ac06b8067c4ac36abfb3d8b2c87bccfe1d9ece35779ca53948a9`, `45ed5c13b4cbab62c511a7300eec057ad2e777d18b6c3400a59b866c440210d2`, `15eeeb0b121b66c7066d64c203ded6f452b192913eb43bb1453fd260985a0252`
- corrected draft aggregate JSON/TSV SHA-256: `a7e7bafa6569f4815aba29f98dc8defa6c7a606e459af2b7913c251b16fca0fb`, `d9d7ad387ecdcc2898f4124a752906a8c631a16fad6f996a0581167a93d7509c`
- backup/apply/recovery helper SHA-256: `a4daf13da8868cbaaffc5013a45b27e52dcf0a0ab3f40c29f71930b6d769145e`, `22a7ae02f71b1f8580f516763116fe668f35492f17beb596b3c31a273eb71c51`, `0426845841338c7c9a9f64b8b5c2fa8aac0e90eca7314d3d86968786f4599777`
- prepared manifest/approval-wrapper SHA-256: `5c900eefa81c5f88ad593d0adf7e51388d99097e04ff6476c06a72f8198783c2`, `a5faed395c0b89407a6017a0b27b72776127047d180aac6bc8c45c261eb693c0`
- v2 aggregate/ready/protected-after SHA-256: `def88e5240f653f98c8d78adf7f8dddcdb05bfda30d08bd0fada7df489969208`, `c26de83abdc9496cd1301470918ec39ecca1cf389ef0ae1c6504da1800d1c431`, `0b6b24e3bb71ec53a109799736432c21a4c3e7e7b9ad9c8e57c7a476462d4b78`
