# Target Rediscovery v0

Laptop-only multi-omic target rediscovery (DepMap CRISPR + RNA-seq) with a scientist-friendly UI.

## Quickstart
\\\powershell
# 1) Activate env (after packages finish installing)
.\.venv\Scripts\Activate.ps1

# 2) Build features and rank
python .\src\build_features.py
python .\src\train_and_rank.py

# 3) Explore UI
streamlit run .\ui\app.py
\\\

## Data (place here)
- \data/raw/CRISPRGeneEffect.csv\
- \data/raw/OmicsExpressionTPMLog1HumanAllGenes.csv\
- \data/raw/Mode.csv\ (model metadata; may be named Model.csv in other releases)

## Outputs
- \data/processed/features_v0.csv\
- \eports/v0_ranked_genes.csv\
- \eports/metrics.txt\
- \eports/top_hits_explanations.csv\ (when computed)

## Methods (short)
Composite score: \S = z_mean_ess + 0.5·z_prop_ess + 0.5·z_lfc_expr\. Weak labels (v0) = top 3% by S for sanity-check classifier; replace with OpenTargets in v0.1.

## License
MIT (see LICENSE)
