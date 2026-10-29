"""Export the human x mouse block of the joint SAMap graph, the kernel of the pseudotime transfer.

Called by rule samap_graph_block (workflow/rules/10c_samap_figures.smk). Reads obs and obsp["connectivities"] from
the joint SAMap object (rows = human nuclei, columns = mouse nuclei, the block samap_transfer_pt.py uses) and writes
C_hs_mm.npz, hs_names.txt, mm_names.txt and meta.json (shapes, nnz, sha256 of the source file) in --out-dir.
"""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp
from anndata.io import read_elem


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--joint", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Only obs and the graph are read; .X (SAM-preprocessed expression) is not needed.
    with h5py.File(args.joint, "r") as f:
        obs = read_elem(f["obs"])
        conn = sp.csr_matrix(read_elem(f["obsp/connectivities"]))
        mapping_k = f["uns/mapping_K"][()].item() if "uns/mapping_K" in f else None

    names = obs.index.astype(str).to_numpy()
    species = obs["species"].astype(str).to_numpy()
    is_human = species == "hs"
    is_mouse = species == "mm"

    C = conn[is_human][:, is_mouse].tocsr()
    sp.save_npz(args.out_dir / "C_hs_mm.npz", C, compressed=True)
    (args.out_dir / "hs_names.txt").write_text("\n".join(names[is_human]) + "\n")
    (args.out_dir / "mm_names.txt").write_text("\n".join(names[is_mouse]) + "\n")

    row_nnz = np.diff(C.indptr)
    meta = {
        "joint": str(args.joint),
        "joint_bytes": args.joint.stat().st_size,
        "joint_sha256": sha256(args.joint),
        "mapping_K": mapping_k,
        "n_human": int(is_human.sum()),
        "n_mouse": int(is_mouse.sum()),
        "nnz": int(C.nnz),
        "row_nnz_min": int(row_nnz.min()),
        "row_nnz_median": float(np.median(row_nnz)),
        "row_nnz_max": int(row_nnz.max()),
        "rows_without_mouse_edge": int((row_nnz == 0).sum()),
    }
    (args.out_dir / "meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
