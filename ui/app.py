import streamlit as st
import pandas as pd
from pathlib import Path

REPORTS = Path("reports")
rank_path = REPORTS / "v0_ranked_genes.csv"
metric_path = REPORTS / "metrics.txt"

st.set_page_config(page_title="Target Rediscovery v0", layout="wide")
st.title("Target Rediscovery v0")

if metric_path.exists():
    st.caption(metric_path.read_text(encoding="utf-8"))
else:
    st.warning("Run the pipeline first to generate metrics.")

if not rank_path.exists():
    st.error("No ranked results found. Run build_features.py and train_and_rank.py, then refresh.")
else:
    df = pd.read_csv(rank_path)
    st.subheader("Top-ranked genes")
    top_k = st.slider("Top-K genes to view", min_value=10, max_value=200, value=50, step=10)
    st.dataframe(df.head(top_k), use_container_width=True)

    st.subheader("Download")
    st.download_button(
        "Download Top-K CSV",
        df.head(top_k).to_csv(index=False).encode("utf-8"),
        "topk.csv",
        "text/csv"
    )
