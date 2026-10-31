"""Draw the Slide-tags sample table of Fig. S12A: stage, section thickness, 10X lanes, nuclei, % mapped.

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for S12A.
Reads cells.csv.gz, lanes.csv and extract_provenance.json from --tables-dir.
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.
"""

import logging

import matplotlib.pyplot as plt
import pandas as pd

import mouse_tissue_style as style


# Half of a 180 mm page, as in the manuscript's S12 layout.
FIG_SIZE_IN = (3.6, 1.9)

# Section thickness is recorded nowhere in the repository (config, sample sheet, docs); these
# are the values of the manuscript's S12A table. Uninduced: two 30 µm sections on one array.
SECTION_THICKNESS = {
    "uninduced": "2 × 30 µm",
    "1.5_dpi_early": "20 µm",
    "1.5_dpi_late": "20 µm",
    "2.5_dpi": "20 µm",
    "3.5_dpi": "20 µm",
    "4.5_dpi_bleeding": "20 µm",
}

# Column header, x position (axes fraction) and alignment, left to right as in the manuscript.
COLUMNS = [
    ("timepoint", 0.0, "left"),
    ("section\nthickness", 0.38, "center"),
    ("#10X\nlanes", 0.545, "center"),
    ("#cells", 0.75, "right"),
    ("spatially\nmapped (%)", 0.915, "center"),
]


def table_rows(cells: pd.DataFrame, lanes: pd.DataFrame) -> pd.DataFrame:
    """One row per stage: nuclei released, 10X lanes and the share with Slide-tags coordinates."""
    per_stage = cells.groupby("timepoint").agg(
        n_nuclei=("has_coordinates", "size"),
        n_mapped=("has_coordinates", "sum"),
    )
    per_stage["pct_mapped"] = 100 * per_stage["n_mapped"] / per_stage["n_nuclei"]
    per_stage = per_stage.join(lanes.set_index("timepoint")[["n_lanes", "samples"]])
    missing = per_stage.index[per_stage["n_lanes"].isna()].tolist()
    if missing:
        raise ValueError(f"stages without sample-sheet lanes: {missing}")
    return per_stage.reindex(style.STAGE_ORDER)


def main() -> None:
    args = style.panel_arguments(__doc__).parse_args()
    style.apply_panel_style()

    cells = style.read_cells(args.tables_dir)
    lanes = pd.read_csv(args.tables_dir / "lanes.csv")
    rows = table_rows(cells, lanes)
    logging.info("Table rows:\n%s", rows.to_string())

    # Figure ----

    # A typeset table: text placed on a blank axes, one line per stage, rules above and below
    # the header as in the manuscript.
    fig, ax = plt.subplots(figsize=FIG_SIZE_IN)
    ax.set_axis_off()
    # Slightly wider than the last column's centre plus half its header, so the rules span it.
    ax.set_xlim(0, 1.03)
    n_rows = len(rows)
    ax.set_ylim(-0.1, n_rows + 2.0)

    # Headers sit on one baseline, so one-line headers align with the second line of the others.
    for header, x, align in COLUMNS:
        ax.text(x, n_rows + 0.3, header, ha=align, va="bottom", fontweight="bold")

    ax.axhline(n_rows + 1.95, color="black", linewidth=1.0)
    ax.axhline(n_rows + 0.12, color="black", linewidth=0.6)

    for i, (stage, row) in enumerate(rows.iterrows()):
        y = n_rows - 0.5 - i
        values = [
            style.TABLE_STAGE_LABELS[stage],
            SECTION_THICKNESS[stage],
            f"{int(row['n_lanes'])}",
            f"{int(row['n_nuclei']):,}",
            f"{row['pct_mapped']:.1f}",
        ]
        for value, (_, x, align) in zip(values, COLUMNS):
            ax.text(x, y, value, ha=align, va="center")

    # Save ----

    provenance = style.base_provenance(
        "S12A",
        args,
        inputs=[args.tables_dir / "cells.csv.gz", args.tables_dir / "lanes.csv"],
        script=__file__,
    )
    provenance.update({
        "definitions": {
            "#cells": "released nuclei per stage (all of them, embedding outliers included)",
            "#10X lanes": "config/samples.tsv rows per stage (mouse, slidetags), one row per lane",
            "spatially mapped (%)": "nuclei with Slide-tags coordinates (x/y_um_dbscan) before spatial QC",
            "section thickness": "manuscript Fig. S12A table; not recorded in the repository",
        },
        "rows": {
            stage: {
                "n_nuclei": int(row["n_nuclei"]),
                "n_mapped": int(row["n_mapped"]),
                "pct_mapped": round(float(row["pct_mapped"]), 2),
                "n_lanes": int(row["n_lanes"]),
                "lanes": row["samples"],
                "section_thickness": SECTION_THICKNESS[stage],
            }
            for stage, row in rows.iterrows()
        },
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
