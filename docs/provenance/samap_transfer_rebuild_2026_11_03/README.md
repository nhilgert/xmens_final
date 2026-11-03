# SAMap pseudotime transfer rebuilt by its rules (3 Nov 2026)

The transfer of record, `code/samap-scvi/full_fibdec/runs_raw/vanilla_raw/pt_transfer_harmonized/` and
`samap_transferred_pt_harmonized.h5ad` (sha256 `62fb75d7…`), was made by hand on 24 Oct 2026 with no log. Its
rules, `samap_raw_pt_transfer` and `samap_raw_attach_pt` (`workflow/rules/21_samap_raw.smk`, variant `_harmonized`),
were run with the rule file, scripts and config of commit `fd1859c` into a scratch root, so the originals were never
touched:

```bash
snakemake <root>/runs_raw/vanilla_raw/samap_transferred_pt_harmonized.h5ad \
    --configfile transfer_overlay.yaml \
    --allowed-rules samap_raw_pt_transfer samap_raw_attach_pt \
    --cores 4 --resources mem_mb=128000
```

`transfer_overlay.yaml` set only `samap_raw.root` to
`tmp/2026_11_03_fig6c_pipeline/samap_raw_root`, which held a link to the joint export
`runs_raw/vanilla_raw/samap_transferred.h5ad`. The mouse pseudotime came from
`data/deliverables/2026_10_22_flagged_harmonized` (column `pt_harmonized`), as the config names it.

| output | result |
|---|---|
| `pt_transferred_hs.csv` | byte-identical (sha256 `5f3535d2…`) |
| `summary.json` | differs only in the two input paths (`joint`, `pt_source`) |
| `samap_transferred_pt_harmonized.h5ad` | every obs column identical, the four transfer columns included (max \|diff\| 0, 47,352 human nuclei with a value); every dataset outside obs identical except `uns/provenance`, whose `pt_transfer.attached` timestamp and `source_csv` path differ |

Files: `transfer_rebuild_report.json` (written by `compare_transfer_rebuild.py`), `transfer_rebuild_snakemake.log`
(the run), `uns_provenance.txt` (both provenance strings). The rebuilt files stay in
`/n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_03_fig6c_pipeline/samap_raw_root/`.
