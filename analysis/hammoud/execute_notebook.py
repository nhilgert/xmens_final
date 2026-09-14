#!/usr/bin/env python
"""
Execute a notebook in-process and write an executed copy with real outputs.

Used because this sandbox forbids binding the local sockets a Jupyter kernel needs
(`PermissionError` on `socket.bind`), so `nbconvert --execute` cannot start one. On a
host where it can, the canonical command is equivalent and preferred:

    python -m jupyter nbconvert --to notebook --execute <nb>.ipynb --output <nb>_executed.ipynb

Text (prints and `display(...)`) is captured as one ordered stream per cell; figures are
captured after the cell body and attached as PNG outputs.
"""
import base64, io, sys, contextlib, traceback
import nbformat as nbf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else "hammoud_figure_walkthrough.ipynb"
DST = sys.argv[2] if len(sys.argv) > 2 else SRC.replace(".ipynb", "_executed.ipynb")

nb = nbf.read(SRC, as_version=4)
ns = {}
plt.show = lambda *a, **k: None          # figures are harvested after each cell

def _display(*objs):
    for o in objs:
        print(o.to_string() if hasattr(o, "to_string") else repr(o))
ns["display"] = _display

n_err = 0
for idx, cell in enumerate(nb.cells):
    if cell.cell_type != "code":
        continue
    buf, outs = io.StringIO(), []
    err = None
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            exec(cell.source, ns)
        except Exception:
            err = traceback.format_exc()
    txt = buf.getvalue()
    if txt.strip():
        outs.append(nbf.v4.new_output("stream", name="stdout", text=txt))
    for num in plt.get_fignums():
        fig = plt.figure(num)
        b = io.BytesIO()
        fig.savefig(b, format="png", dpi=110, bbox_inches="tight")
        outs.append(nbf.v4.new_output(
            "display_data",
            data={"image/png": base64.b64encode(b.getvalue()).decode(),
                  "text/plain": f"<Figure {fig.get_size_inches()[0]:.1f}x{fig.get_size_inches()[1]:.1f}>"},
            metadata={}))
    plt.close("all")
    if err:
        n_err += 1
        outs.append(nbf.v4.new_output("error", ename="Error", evalue="see traceback",
                                      traceback=err.splitlines()))
        print(f"[cell {idx}] ERROR\n{err[-1400:]}", file=sys.stderr)
    cell.outputs = outs
    cell.execution_count = idx

nbf.write(nb, DST)
n_fig = sum(1 for c in nb.cells for o in (c.get("outputs") or [])
            if "image/png" in (o.get("data") or {}))
print(f"wrote {DST}: {len(nb.cells)} cells, {n_fig} figures, {n_err} cell errors")
sys.exit(1 if n_err else 0)
