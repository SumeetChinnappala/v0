import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression

PROC = Path("data/processed")
REPORTS = Path("reports"); REPORTS.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(PROC/"features_v0.csv")  # columns: gene, z_mean_ess, z_prop_ess, z_lfc_expr, S

# Weak labels for v0: use top 3% by S as positives (replace with OpenTargets later)
k = max(1, int(0.03*len(df)))
y = np.zeros(len(df), dtype=int)
y[df["S"].nlargest(k).index] = 1

X = df[["z_mean_ess","z_prop_ess","z_lfc_expr"]].values
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
clf = LogisticRegression(max_iter=2000)

probs = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:,1]
auroc = roc_auc_score(y, probs)
auprc = average_precision_score(y, probs)

out = df.copy()
out["p_known_like"] = probs
out = out.sort_values(["S","p_known_like"], ascending=False)
out.to_csv(REPORTS/"v0_ranked_genes.csv", index=False)

with open(REPORTS/"metrics.txt","w") as f:
    f.write(f"AUROC={auroc:.3f}\nAUPRC={auprc:.3f}\n")

print(f"[OK] Saved reports/v0_ranked_genes.csv and metrics (AUROC={auroc:.3f}, AUPRC={auprc:.3f})")
