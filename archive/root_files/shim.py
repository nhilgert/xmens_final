
import yaml, types, pathlib, builtins
X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
cfg = yaml.safe_load(open(X/"config/config.yaml"))
class NS(dict): __getattr__ = dict.__getitem__
snakemake = types.SimpleNamespace(
    input=NS(h5ad=str(X/"data/interim/annotated_mouse.h5ad")),
    output=NS(h5ad=str(X/"data/interim/qc_spatial_mouse.h5ad"), tally=str(X/"results/tables/spatial_qc_mouse.csv")),
    params=NS(sqc=cfg["spatial_qc"]), wildcards=NS(species="mouse"), log=[str(X/"logs/04f_spatial_qc_mouse.log")])
builtins.snakemake = snakemake
exec(compile(open(X/"workflow/scripts/spatial_qc.py").read(), str(X/"workflow/scripts/spatial_qc.py"), "exec"), {"snakemake": snakemake, "__name__": "__main__"})
