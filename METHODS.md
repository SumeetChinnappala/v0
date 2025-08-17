# METHODS (v0)

**Data release:** DepMap Public (release + date), files: CRISPRGeneEffect.csv, OmicsExpressionTPMLog1HumanAllGenes.csv, Mode/Model.csv.

**Preprocess:**
- Align by DepMap_ID intersection between CRISPR and RNA-seq.
- Leukemia lineage: \Blood\, \Haematopoietic_and_lymphoid_tissue\.

**Features:**
- \z_mean_ess\: invert CRISPR gene effect (lower = more essential).
- \z_prop_ess\: fraction of leukemia lines with gene effect < -0.5.
- \z_lfc_expr\: mean(leukemia) - mean(non-leukemia) on log(TPM+1).

**Score:** \S = z_mean_ess + 0.5*z_prop_ess + 0.5*z_lfc_expr\.

**Validation (v0):**
- Precision/Recall@K for canonical leukemia genes (ABL1, FLT3, BCL2, JAK2, KIT, IDH1/2, DNMT3A).
- Enrichment (Fisher's exact, top decile).
- Ablation: multi-omic vs CRISPR-only vs RNA-only.

**Runtime & hardware:** Windows 11 laptop, CPU-only.
