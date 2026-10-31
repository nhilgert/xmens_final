"""Fig. S17Aatlas: the fibroblast/decidual nuclei of Fig. S17A at their whole-tissue UMAP coordinates.

Reads the fib/dec plot table (subtypes) and the whole-tissue UMAP table (obsm X_umap_scVI of the release's
whole-tissue object, with its embedding-outlier flags), both written by extract_panel_table.py, joins them on
obs_names and draws each fib/dec nucleus where it sits in the whole-tissue embedding, coloured by subtype with
the S17A colours and direct labels. Writes <out>.svg (live text), <out>.png (300 dpi) and <out>.provenance.json.
Rule: panels_mouse_stroma_s17aatlas (workflow/rules/91_panels_mouse_stroma.smk).
"""
import argparse
import datetime
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import panel_style as ps

SEED = 0  # draw order is shuffled as in S17A, so no subtype systematically covers another
SVG_HASH_SALT = "panel_S17Aatlas"  # fixed SVG element ids: the same tables give the same bytes

# S17A's width and marker area; the height follows the window's aspect (equal axis scales). The embedding is
# drawn as stored, as in Fig. 3A (the same whole-tissue fit), so the two panels share their orientation.
FIG_WIDTH_IN = 2.7
MARKER_AREA = 0.6  # pt^2

# Window around the drawn nuclei (UMAP units): room for the labels set outside the band.
WINDOW_PAD = {"left": 0.3, "right": 0.3, "bottom": 0.3, "top": 1.2}

# Labels: text position relative to the subtype's median (UMAP units), horizontal alignment, and the point a
# leader runs to (relative to the median), or None for a label beside its cloud. Chosen by eye on the
# rendered panel: names sit in the white space next to their own cloud; eFib2i is enclosed by eFib2, eFib1,
# eFibP and dec1, so its name sits above the band with a leader, as dec1's does in S17A.
ATLAS_LABELS = {
    "eFib3": ((-0.31, 2.72), "center", None),
    "eFib2": ((0.30, 2.31), "center", None),
    "eFib2i": ((0.81, 3.54), "center", (0.81, 1.69)),
    "dec1": ((0.34, 2.11), "center", None),
    "dec2": ((-0.76, -3.92), "center", None),
    "dec3": ((0.98, -3.44), "center", None),
    "eFib1": ((-1.37, -0.22), "right", None),
    "eFibP": ((-0.32, -1.84), "center", None),
    "decAp": ((2.11, -0.09), "left", None),
}


def joined_nuclei(fibdec: pd.DataFrame, atlas: pd.DataFrame) -> pd.DataFrame:
    """Every fib/dec nucleus with its whole-tissue coordinates and both outlier flags, matched by obs_names."""
    missing = fibdec.index.difference(atlas.index)
    if len(missing):
        raise KeyError(f"{len(missing)} fib/dec nuclei are not in the whole-tissue table, e.g. {list(missing[:3])}")
    joined = fibdec[["cell_subtype", "embedding_outlier"]].rename(columns={"embedding_outlier": "outlier_fibdec"})
    joined = joined.join(atlas[["X_umap_scVI_1", "X_umap_scVI_2", "cell_type", "embedding_outlier"]])
    return joined.rename(columns={"embedding_outlier": "outlier_atlas", "X_umap_scVI_1": "x", "X_umap_scVI_2": "y"})


def save(fig, stem: Path, provenance: dict) -> None:
    """Write <stem>.svg and .png byte-reproducibly (no dates, fixed SVG ids), then the provenance."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    provenance["text_overlaps"] = ps.text_overlaps(fig)
    files = {}
    for ext, metadata in (("svg", {"Date": None}), ("png", None)):
        path = stem.parent / f"{stem.name}.{ext}"
        # The scatter is rasterized at the PNG resolution inside the SVG; text stays live.
        fig.savefig(path, dpi=ps.PNG_DPI, metadata=metadata)
        files[path.name] = ps.sha256_of(path)
    plt.close(fig)
    provenance["written"] = datetime.datetime.now().isoformat(timespec="seconds")
    provenance["matplotlib"] = matplotlib.__version__
    provenance["font_family"] = ps.FONT_FAMILY
    provenance["files"] = files
    (stem.parent / f"{stem.name}.provenance.json").write_text(json.dumps(provenance, indent=1, default=str))
    print(f"wrote {stem}.{{svg,png,provenance.json}}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path, help="fib/dec plot table")
    parser.add_argument("--table-prov", required=True, type=Path)
    parser.add_argument("--atlas-table", required=True, type=Path, help="whole-tissue UMAP table")
    parser.add_argument("--atlas-table-prov", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="output stem")
    args = parser.parse_args()

    # Load ----

    fibdec, fibdec_record = ps.load_table(args.table, args.table_prov)
    atlas, atlas_record = ps.load_table(args.atlas_table, args.atlas_table_prov)
    unknown = set(fibdec["cell_subtype"]) - set(ps.SUBTYPES)
    if unknown:
        raise ValueError(f"subtypes outside the S17 vocabulary: {sorted(unknown)}")
    nuclei = joined_nuclei(fibdec, atlas)

    # The whole-tissue object's flags belong to the embedding drawn here; the release copies them into the
    # fib/dec object, so the two should agree (counted in the provenance).
    shown = nuclei[~nuclei["outlier_atlas"]]
    print(f"{len(nuclei)} fib/dec nuclei matched; {int(nuclei['outlier_atlas'].sum())} embedding outliers hidden")
    order = np.random.default_rng(SEED).permutation(len(shown))

    # Figure ----

    ps.apply_style(dark=False)
    plt.rcParams["svg.hashsalt"] = SVG_HASH_SALT
    x_lim = (shown["x"].min() - WINDOW_PAD["left"], shown["x"].max() + WINDOW_PAD["right"])
    y_lim = (shown["y"].min() - WINDOW_PAD["bottom"], shown["y"].max() + WINDOW_PAD["top"])
    height_in = FIG_WIDTH_IN * (y_lim[1] - y_lim[0]) / (x_lim[1] - x_lim[0])
    fig = plt.figure(figsize=(FIG_WIDTH_IN, height_in))
    ax = fig.add_axes([0, 0, 1, 1])

    # One point per nucleus at its whole-tissue position, colour = subtype as in S17A; shuffled draw order.
    ax.scatter(
        shown["x"].to_numpy()[order],
        shown["y"].to_numpy()[order],
        c=shown["cell_subtype"].map(ps.SUBTYPE_COLORS).to_numpy()[order],
        s=MARKER_AREA,
        linewidths=0,
        rasterized=True,
    )

    # Direct labels in the S17A style: black, title size, beside each cloud; labels of clouds enclosed by
    # others sit outside the band with a leader to their cloud.
    medians = shown.groupby("cell_subtype")[["x", "y"]].median()
    label_record = {}
    for subtype in ps.SUBTYPES:
        (dx, dy), ha, leader_to = ATLAS_LABELS[subtype]
        anchor = medians.loc[subtype]
        position = (anchor["x"] + dx, anchor["y"] + dy)
        if leader_to is None:
            ax.text(*position, subtype, ha=ha, va="center", fontsize=ps.SIZE_TITLE)
        else:
            target = (anchor["x"] + leader_to[0], anchor["y"] + leader_to[1])
            ax.annotate(subtype, xy=target, xytext=position, ha=ha, va="center", fontsize=ps.SIZE_TITLE,
                        arrowprops=dict(arrowstyle="-", lw=0.8, color="black", shrinkA=1, shrinkB=0))
        label_record[subtype] = {"text_xy": [round(v, 3) for v in position], "leader": leader_to is not None}

    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_aspect("equal")
    ax.axis("off")

    # Save ----

    agree = int((nuclei["outlier_atlas"] == nuclei["outlier_fibdec"]).sum())
    provenance = {
        "panel": "S17Aatlas",
        "script": Path(__file__).name,
        "inputs": {"fibdec": fibdec_record, "whole_tissue": atlas_record},
        "embedding": "obsm['X_umap_scVI'] of the whole-tissue object (whole-tissue fit), as stored, as in Fig. 3A",
        "matching": f"{len(nuclei)} of {len(fibdec)} fib/dec nuclei found in the whole-tissue table by obs_names",
        "whole_tissue_cell_type_of_fibdec_nuclei": nuclei["cell_type"].value_counts().to_dict(),
        "embedding_outliers_hidden": int(nuclei["outlier_atlas"].sum()),
        "outlier_flag_source": "whole-tissue object (the embedding drawn)",
        "outlier_flags_agree_with_fibdec_object": f"{agree} of {len(nuclei)}",
        "nuclei_drawn": int(len(shown)),
        "nuclei_by_subtype_drawn": shown["cell_subtype"].value_counts().reindex(ps.SUBTYPES).to_dict(),
        "window_umap_units": {"x": [round(v, 3) for v in x_lim], "y": [round(v, 3) for v in y_lim]},
        "labels": label_record,
        "colours": ps.SUBTYPE_COLORS,
        "colour_source": "as Fig. S17A (sampled from the manuscript Fig. S17B bars)",
        "draw_order": f"numpy default_rng({SEED}) permutation",
        "pending_labels": ps.pending_label_notes(shown["cell_subtype"]),
    }
    save(fig, args.out, provenance)


if __name__ == "__main__":
    main()
