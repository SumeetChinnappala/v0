import os, re, pandas as pd, numpy as np
from pathlib import Path

RAW = Path("data/raw")
PROC = Path("data/processed"); PROC.mkdir(parents=True, exist_ok=True)

def find_one(patterns):
    for p in RAW.glob("*.csv"):
        name = p.name.lower()
        if any(all(tok in name for tok in pat) for pat in patterns):
            return p
    return None

f_ge   = find_one([["crispr","gene","effect"]])
f_expr = find_one([["tpm","log1"],["tpm","logp1"],["tpmlog1"]])
f_meta = find_one([["model"],["mode"]])
if not f_ge or not f_expr or not f_meta:
    raise FileNotFoundError(f"Missing files: CRISPR={f_ge}, TPM={f_expr}, META={f_meta}")

ge   = pd.read_csv(f_ge)
expr = pd.read_csv(f_expr)
meta = pd.read_csv(f_meta)

def normalize_gene_df(df):
    for c in ["gene","Gene","gene_symbol","Hugo_Symbol","HUGO_symbol","Unnamed: 0","symbol"]:
        if c in df.columns:
            out = df.set_index(c)
            return out[~out.index.duplicated(keep="first")]
    first = df.columns[0]
    if re.search(r"depmap|model", first, flags=re.I):
        raise ValueError("No gene column found in matrix")
    out = df.set_index(first)
    return out[~out.index.duplicated(keep="first")]

ge_i   = normalize_gene_df(ge)
expr_i = normalize_gene_df(expr)

common_cells = sorted(set(ge_i.columns).intersection(expr_i.columns))
if not common_cells:
    raise RuntimeError("No overlapping cell lines between CRISPR and expression.")

meta_cols = {c.lower(): c for c in meta.columns}
def pick(*opts):
    for o in opts:
        if o in meta_cols: return meta_cols[o]
    for o in opts:
        if o in meta.columns: return o
    return None

cid  = pick("depmap_id","model_id","modelid","model","model_id_broad","ModelID","Model")
clin = pick("lineage","lineage_mapped","primary_lineage","Lineage","primary_disease")
if not cid:  raise RuntimeError("Model metadata missing DepMap_ID/ModelID column")
if not clin: raise RuntimeError("Model metadata missing lineage column")

meta = meta.rename(columns={cid:"DepMap_ID", clin:"lineage"})
meta = meta[meta["DepMap_ID"].isin(common_cells)].copy()
meta["lineage"] = meta["lineage"].astype(str)

LEUK = {"Blood","Haematopoietic_and_lymphoid_tissue","Haematopoietic and lymphoid tissue"}
leuk = meta.loc[meta["lineage"].isin(LEUK), "DepMap_ID"].tolist()
nonl = [c for c in common_cells if c not in set(leuk)]

ge_c, expr_c = ge_i[common_cells], expr_i[common_cells]

g_mean_ess = ge_c[leuk].mean(1) if leuk else ge_c.mean(1)
g_prop_ess = (ge_c[leuk] < -0.5).mean(1) if leuk else (ge_c < -0.5).mean(1)
lfc_expr   = (expr_c[leuk].mean(1) - (expr_c[nonl].mean(1) if nonl else expr_c[leuk].mean(1)))

def z(s):
    s = s.astype(float); m = s.mean(); sd = s.std(ddof=0)
    return pd.Series(0.0, index=s.index) if (sd==0 or not np.isfinite(sd)) else (s-m)/sd

feat = pd.DataFrame({
  "z_mean_ess": z(-g_mean_ess),
  "z_prop_ess": z(g_prop_ess),
  "z_lfc_expr": z(lfc_expr),
}).fillna(0.0)
feat["S"] = feat.eval("z_mean_ess + 0.5*z_prop_ess + 0.5*z_lfc_expr")

feat.index.name = "gene"
feat.reset_index().to_csv(PROC/"features_v0.csv", index=False)
print("[OK] Wrote data/processed/features_v0.csv")
