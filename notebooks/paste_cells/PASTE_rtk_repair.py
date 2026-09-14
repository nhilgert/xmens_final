# ==========================================================================
# PASTE THIS over the failing RealTimeKernel cell, then run it.
# Your `tp` (moscot solve) is already in the kernel — this does NOT re-solve.
#
# Why: cellrank's transition_matrix setter demands every row sum to 1 within
# rtol=1e-12. The stitched moscot couplings leave a few all-zero rows (cells no
# coupling reached — the "contains 1 empty row(s)" warnings you saw at (1.5,2.0)
# and (3.0,4.0)). This is the workaround from cells 28-30 of your original
# realtime_informed_pseudotime_human.ipynb, which is load-bearing, not diagnostic.
# ==========================================================================
import numpy as np
from scipy.sparse import csr_matrix

rtk = cr.kernels.RealTimeKernel.from_moscot(tp)

try:
    rtk = rtk.compute_transition_matrix(
        self_transitions=CONFIG["rtk_self_transitions"],
        conn_weight=CONFIG["rtk_conn_weight"],
        threshold=CONFIG["rtk_threshold"],
    )
    T = rtk.transition_matrix
    print("RealTimeKernel accepted the stitched matrix directly")

except ValueError as e:
    if "row stochastic" not in str(e):
        raise
    print(f"RealTimeKernel rejected the matrix: {e}")
    print("applying the documented repair ...")

    # grab the matrix the setter refused, or re-stitch without going through it
    if getattr(rtk, "_transition_matrix", None) is not None:
        T = rtk._transition_matrix
    else:
        tmap = rtk._restich_couplings(
            rtk.couplings,
            self_transitions=CONFIG["rtk_self_transitions"],
            conn_weight=CONFIG["rtk_conn_weight"],
        )
        T = tmap.X

    T = csr_matrix(T).tocsr()
    rs = np.asarray(T.sum(axis=1)).ravel()
    dead = np.where(~np.isfinite(rs) | (rs == 0))[0]
    print(f"  dead rows (zero / non-finite): {len(dead)}")

    # a dead row becomes uniform over cells at the SAME timepoint
    day_vals = adata.obs["day"].astype(float).to_numpy()
    if len(dead):
        T_lil = T.tolil()
        for i in dead:
            same_tp = np.where(day_vals == day_vals[i])[0]
            T_lil.rows[i] = list(same_tp)
            T_lil.data[i] = [1.0 / len(same_tp)] * len(same_tp)
        T = T_lil.tocsr()

    # renormalize EVERY row to exactly 1.0
    rs = np.asarray(T.sum(axis=1)).ravel()
    T = csr_matrix(T.multiply(1.0 / rs[:, None]))
    rs2 = np.asarray(T.sum(axis=1)).ravel()
    print(f"  after repair: row sums [{rs2.min():.15f}, {rs2.max():.15f}]")
    print(f"  non-stochastic (rtol=1e-12): {int((~np.isclose(rs2, 1.0, rtol=1e-12)).sum())}")

    rtk.transition_matrix = T
    T = rtk.transition_matrix
    print("  assigned successfully")

rs = np.asarray(T.sum(axis=1)).ravel()
print(f"\ntransition matrix {T.shape}  row sums [{rs.min():.9f}, {rs.max():.9f}]")
assert np.allclose(rs, 1.0, atol=1e-6)
