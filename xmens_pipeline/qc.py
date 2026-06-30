"""Stage 01 — quality control, per Cevrim et al. 2025 Methods.

Input: per-sample AnnData from 00b (X = CellBender-corrected counts;
obs has species/timepoint/coords; coords nullable; all cells present).

Procedure (Scanpy):
  1. Compute QC metrics (mito fraction etc.).
  2. Minimum 500 corrected UMIs/cell.
  3. Adaptive 3-MAD filters on log1p(total_counts) and log1p(n_genes_by_counts).
  4. Mitochondrial fraction > 8% excluded.
  5. DoubletDetection (BoostClassifier, 40 iters, voter_thresh 0.5); remove doublets.
     Requires the doubletdetection package (separate env via 01_qc.smk `conda:`);
     if unavailable, doublet removal is SKIPPED with a warning.

Writes data/interim/qc/{sample}.h5ad. Corrected counts kept in layers['counts'].
"""
import numpy as np
import scanpy as sc

inp = snakemake.input.h5ad           # noqa: F821 (snakemake-injected)
out = snakemake.output.h5ad           # noqa: F821
sample = snakemake.wildcards.sample   # noqa: F821
species = snakemake.params.species    # noqa: F821
cfg = snakemake.config                # noqa: F821
logfile = snakemake.log[0]            # noqa: F821

qc = cfg["qc"]
seed = cfg.get("seed", 0)
mito_prefix = cfg["mito_prefix"][species]
dd_cfg = cfg.get("doublet", {})


def mad_bounds(x, n_mads):
    """Lower/upper bound at n_mads scaled-MADs from the median."""
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    scaled = 1.4826 * mad   # consistent estimator of sigma under normality
    return med - n_mads * scaled, med + n_mads * scaled


with open(logfile, "w") as log:
    def say(msg):
        print(f"[qc:{sample}] {msg}", file=log, flush=True)

    np.random.seed(seed)

    adata = sc.read_h5ad(inp)
    n0 = adata.n_obs
    say(f"loaded {n0} cells x {adata.n_vars} genes (species={species})")

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # --- 1. QC metrics ---
    adata.var["mito"] = adata.var_names.str.startswith(mito_prefix)
    n_mito = int(adata.var["mito"].sum())
    say(f"flagged {n_mito} mito genes (prefix '{mito_prefix}')")
    if n_mito == 0:
        say("WARNING: no mito genes matched - check prefix/gene naming for this species")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito"], percent_top=None, log1p=True, inplace=True
    )

    # --- 2. min UMI floor ---
    min_counts = qc["min_counts"]
    keep_umi = adata.obs["total_counts"] >= min_counts
    say(f"min {min_counts} UMIs: keep {int(keep_umi.sum())}/{n0}")

    # --- 3. adaptive 3-MAD on log1p(total_counts) and log1p(n_genes) ---
    n_mads = qc["n_mads"]
    log_counts = np.log1p(adata.obs["total_counts"].to_numpy())
    log_genes = np.log1p(adata.obs["n_genes_by_counts"].to_numpy())
    c_lo, c_hi = mad_bounds(log_counts, n_mads)
    g_lo, g_hi = mad_bounds(log_genes, n_mads)
    keep_counts = (log_counts >= c_lo) & (log_counts <= c_hi)
    keep_genes = (log_genes >= g_lo) & (log_genes <= g_hi)
    say(f"{n_mads}-MAD log-counts keep {int(keep_counts.sum())}; "
        f"log-genes keep {int(keep_genes.sum())}")

    # --- 4. mito fraction ---
    max_mito = qc["max_pct_mito"]
    keep_mito = adata.obs["pct_counts_mito"] <= max_mito
    say(f"mito <= {max_mito}%: keep {int(keep_mito.sum())}")

    keep = keep_umi & keep_counts & keep_genes & keep_mito
    adata = adata[keep.to_numpy()].copy()
    say(f"after threshold QC: {adata.n_obs}/{n0} cells")

    # --- 5. DoubletDetection (optional env) ---
    try:
        import doubletdetection
        ver = getattr(doubletdetection, "__version__", "?")
        say(f"DoubletDetection {ver}: {dd_cfg.get('n_iters', 40)} iters, "
            f"voter_thresh={dd_cfg.get('voter_thresh', 0.5)}")
        clf = doubletdetection.BoostClassifier(
            n_iters=dd_cfg.get("n_iters", 40),
            clustering_algorithm="leiden",
            standard_scaling=True,
            pseudocount=0.1,
            n_jobs=-1,
            random_state=seed,
        )
        counts = adata.layers["counts"]
        labels = clf.fit(counts).predict(
            p_thresh=dd_cfg.get("p_thresh", 1e-7),
            voter_thresh=dd_cfg.get("voter_thresh", 0.5),
        )
        adata.obs["doublet"] = labels.astype(bool)
        adata.obs["doublet_score"] = clf.doublet_score()
        n_doub = int(adata.obs["doublet"].sum())
        say(f"flagged {n_doub} doublets; removing")
        adata = adata[~adata.obs["doublet"].to_numpy()].copy()
    except ImportError:
        say("WARNING: doubletdetection not available in this env - SKIPPING doublet "
            "removal. Set config['envs']['doubletdetection'] to an env with it "
            "(==4.2 to match the paper).")
    except Exception as e:  # noqa: BLE001
        say(f"WARNING: DoubletDetection failed ({e}); proceeding without doublet removal")

    say(f"final: {adata.n_obs} cells x {adata.n_vars} genes "
        f"({n0 - adata.n_obs} removed total)")
    adata.write_h5ad(out)
    say(f"wrote -> {out}")
