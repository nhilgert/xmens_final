"""Reading copy of one species' Leiden review dot plots: one page per annotation step, in review order.

Called by rule panels_marker_dotplots_review_pdf (workflow/rules/91_panels_marker_dotplots.smk), once per
species. Reads the review tables of every step (panels_marker_dotplots_review_extract), Table M1 and the
annotation vocabulary, and draws each page with the code of plot_review_step.py, so the pages match the step
SVGs. Writes one PDF without creation date or document id: a re-render is byte-identical.
"""

import argparse
import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import yaml
from matplotlib.backends.backend_pdf import PdfPages

import dotplot_style as style
import marker_panels as panels
from plot_review_step import REVIEW_RC, draw_step, load_step


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=list(panels.REVIEW_PLOT))
    parser.add_argument("--steps", required=True, nargs="+", help="annotation steps in page order")
    parser.add_argument("--tables-dir", required=True, type=Path)
    parser.add_argument("--table-m1", required=True, type=Path)
    parser.add_argument("--vocabulary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    table_m1 = panels.read_table_m1(args.table_m1)
    vocabulary = yaml.safe_load(args.vocabulary.read_text())
    args.out.parent.mkdir(parents=True, exist_ok=True)

    with mpl.rc_context(REVIEW_RC), PdfPages(args.out, metadata=style.PDF_METADATA) as pdf:
        for step in args.steps:
            data = load_step(args.tables_dir, args.species, step, table_m1, vocabulary)
            fig, _, _, _ = draw_step(args.species, step, data, vocabulary)
            pdf.savefig(fig)
            plt.close(fig)
            logging.info("%s %s: page %d", args.species, step, pdf.get_pagecount())


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    # Font subsetting logs every glyph list at INFO.
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    main()
