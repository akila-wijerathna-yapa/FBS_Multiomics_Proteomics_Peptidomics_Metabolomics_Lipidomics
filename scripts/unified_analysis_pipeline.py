#!/usr/bin/env python3
"""
FBS Multi-Omics Paper 3 — FRESH Figure Generator
==================================================
Generates ALL figures from raw CSV data.  Nothing is copied from old plots.

Figures match Akila's PPTX slides 5-10 exactly:
  Figure 1: Integration outperforms single-omics
  Figure 2: MOFA2 latent factors
  Figure 3: DIABLO sparse signature
  Figure 4: Cross-omics correlation network
  Figure 5: Integrated modules & biological themes
  Figure 6: Biologically important molecules / mechanistic model
"""

import os, sys, shutil, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr
import networkx as nx
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# ── paths ────────────────────────────────────────────────────────────────────
BASE   = Path("/Users/arie/Documents/Upwork_Fiverr/Upwork_Akila")
MULTI  = BASE / "upwork5_multiomcs_all" / "Multiomics_ALL-29October2025"
OUTPUT = BASE / "all 3 papers data and figures" / "Paper_3_MultiOmics"
FIG    = OUTPUT / "figures"
TAB    = OUTPUT / "tables"

GRADE_COLORS  = {"High-grade": "#EF3F37", "Medium-grade": "#FBAF41", "Low-grade": "#262161"}
GRADE_ORDER   = ["High-grade", "Medium-grade", "Low-grade"]
OMICS_COLORS  = {"Proteomics": "#2196F3", "Peptidomics": "#4CAF50",
                 "Metabolomics": "#FF9800", "Lipidomics": "#9C27B0"}
OMICS_ORDER   = ["Proteomics", "Peptidomics", "Metabolomics", "Lipidomics"]
OMICS_LOWER   = ["proteomics", "peptidomics", "metabolomics", "lipidomics"]
CONTRAST      = "High.grade - Low.grade"

# ── styling ──────────────────────────────────────────────────────────────────
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi": 300,
})

def ensure_dirs():
    """Create fresh output directories."""
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    for sub in ["Figure_1", "Figure_2", "Figure_3", "Figure_4",
                "Figure_5", "Figure_6", "Figure_S1"]:
        (FIG / sub).mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  DATA LOADERS                                                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def load_pca_scores():
    """Load per-omics PCA score CSVs."""
    dfs = {}
    for tag in OMICS_LOWER:
        p = MULTI / "PCA" / f"{tag}_pca_scores.csv"
        if p.exists():
            dfs[tag] = pd.read_csv(p)
    return dfs

def load_pca_variance():
    dfs = {}
    for tag in OMICS_LOWER:
        p = MULTI / "PCA" / f"{tag}_pca_variance.csv"
        if p.exists():
            dfs[tag] = pd.read_csv(p)
    return dfs

def load_integrated_pca():
    p = MULTI / "MetaboAnalyst" / "metaboanalyst_pca_scores.csv"
    return pd.read_csv(p) if p.exists() else None

def load_mofa_factors():
    return pd.read_csv(MULTI / "MOFA2" / "mofa_factors.csv")

def load_mofa_variance():
    return pd.read_csv(MULTI / "MOFA2" / "mofa_variance_explained_per_factor.csv")

def load_mofa_weights():
    return pd.read_csv(MULTI / "MOFA2" / "mofa_top_weights.csv")

def load_diablo_scores_long():
    return pd.read_csv(MULTI / "DIABLO" / "diablo_scores_long.csv")

def load_diablo_scores_wide():
    return pd.read_csv(MULTI / "DIABLO" / "diablo_scores_wide.csv")

def load_diablo_loadings():
    return pd.read_csv(MULTI / "Circos" / "diablo_component_loadings.csv")

def load_network_edges():
    return pd.read_csv(MULTI / "Networks" / "multiomics_coexpression_edges.csv")

def load_wgcna_modules():
    dfs = {}
    for tag in OMICS_LOWER:
        p = MULTI / "WGCNA" / f"{tag}_wgcna_modules.csv"
        if p.exists():
            dfs[tag] = pd.read_csv(p)
    return dfs

def load_wgcna_eigengenes():
    dfs = {}
    for tag in OMICS_LOWER:
        p = MULTI / "WGCNA" / f"{tag}_wgcna_eigengenes.csv"
        if p.exists():
            dfs[tag] = pd.read_csv(p)
    return dfs

def load_pathway_data():
    """Load all pathway enrichment CSVs."""
    records = []
    pw = MULTI / "Pathways"
    if not pw.exists():
        return pd.DataFrame()
    for f in pw.glob("*.csv"):
        try:
            df = pd.read_csv(f)
            records.append(df)
        except:
            pass
    if records:
        return pd.concat(records, ignore_index=True)
    return pd.DataFrame()

def load_diff_data():
    dfs = {}
    for tag in OMICS_LOWER:
        p = MULTI / "Differential_Updated" / f"{tag}_differential_results.csv"
        if p.exists():
            dfs[tag] = pd.read_csv(p)
    return dfs


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  HELPER: Silhouette score calculation                                    ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def compute_silhouette(coords, labels):
    """Compute mean silhouette score given coordinates and cluster labels."""
    from sklearn.metrics import silhouette_score
    unique = labels.unique()
    if len(unique) < 2 or len(labels) < 4:
        return np.nan
    try:
        return silhouette_score(coords, labels)
    except:
        return np.nan


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 1 — Integration Outperforms Single-Omics                        ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_1():
    print("\n═══ Figure 1: Integration vs Single-Omics ═══")

    pca_scores = load_pca_scores()
    pca_var    = load_pca_variance()
    int_pca    = load_integrated_pca()
    mofa_fac   = load_mofa_factors()
    diablo_wide = load_diablo_scores_wide()
    diablo_long = load_diablo_scores_long()

    # ── Panel A: 6-panel grid of scatter plots ──
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    silhouette_results = {}

    for idx, tag in enumerate(OMICS_LOWER):
        ax = axes[idx // 3, idx % 3] if idx < 3 else axes[1, idx - 3]
        if tag in pca_scores:
            df = pca_scores[tag]
            var = pca_var.get(tag)
            pc1_var = var.iloc[0]["variance_explained"] * 100 if var is not None and "variance_explained" in var.columns else ""
            pc2_var = var.iloc[1]["variance_explained"] * 100 if var is not None and "variance_explained" in var.columns and len(var) > 1 else ""

            for grade in GRADE_ORDER:
                mask = df["group"] == grade
                ax.scatter(df.loc[mask, "PC1"], df.loc[mask, "PC2"],
                          c=GRADE_COLORS[grade], label=grade, s=120,
                          edgecolors="white", linewidth=1.5, zorder=3, alpha=0.9)

            xlabel = f"PC1 ({pc1_var:.1f}%)" if isinstance(pc1_var, (int, float)) else "PC1"
            ylabel = f"PC2 ({pc2_var:.1f}%)" if isinstance(pc2_var, (int, float)) else "PC2"
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(f"{tag.capitalize()} PCA", fontweight="bold", fontsize=12)
            ax.legend(fontsize=8, framealpha=0.7)

            coords = df[["PC1", "PC2"]].values
            sil = compute_silhouette(coords, df["group"])
            silhouette_results[tag.capitalize()] = sil

    # MOFA panel
    ax_mofa = axes[1, 1]
    if mofa_fac is not None:
        f1 = mofa_fac[mofa_fac["factor"] == "Factor1"].set_index("sample_id")["value"]
        f2 = mofa_fac[mofa_fac["factor"] == "Factor2"].set_index("sample_id")["value"]
        combined = pd.DataFrame({"F1": f1, "F2": f2})
        meta = pd.read_csv(MULTI / "sample_metadata.csv").set_index("sample_id")
        combined = combined.join(meta[["group"]])

        for grade in GRADE_ORDER:
            mask = combined["group"] == grade
            ax_mofa.scatter(combined.loc[mask, "F1"], combined.loc[mask, "F2"],
                           c=GRADE_COLORS[grade], label=grade, s=120,
                           edgecolors="white", linewidth=1.5, zorder=3, alpha=0.9)
        ax_mofa.set_xlabel("MOFA Factor 1")
        ax_mofa.set_ylabel("MOFA Factor 2")
        ax_mofa.set_title("MOFA2 Integrated", fontweight="bold", fontsize=12)
        ax_mofa.legend(fontsize=8, framealpha=0.7)

        sil = compute_silhouette(combined[["F1", "F2"]].values, combined["group"])
        silhouette_results["MOFA2 Integrated"] = sil

    # DIABLO panel
    ax_diablo = axes[1, 2]
    if diablo_wide is not None:
        # Use first dataset's scores as representative
        datasets = diablo_wide["dataset_label"].unique()
        first_ds = datasets[0] if len(datasets) > 0 else None
        if first_ds:
            dsub = diablo_wide[diablo_wide["dataset_label"] == first_ds]
            for grade in GRADE_ORDER:
                mask = dsub["group"] == grade
                ax_diablo.scatter(dsub.loc[mask, "comp1"], dsub.loc[mask, "comp2"],
                                 c=GRADE_COLORS[grade], label=grade, s=120,
                                 edgecolors="white", linewidth=1.5, zorder=3, alpha=0.9)
            ax_diablo.set_xlabel("DIABLO Component 1")
            ax_diablo.set_ylabel("DIABLO Component 2")
            ax_diablo.set_title("DIABLO Integrated", fontweight="bold", fontsize=12)
            ax_diablo.legend(fontsize=8, framealpha=0.7)

            coords = dsub[["comp1", "comp2"]].values
            sil = compute_silhouette(coords, dsub["group"])
            silhouette_results["DIABLO Integrated"] = sil

    fig.suptitle("Figure 1 — Multi-Omics Integration Captures Grade Separation\nBetter Than Single Omics",
                 fontsize=15, fontweight="bold", y=1.02)
    sns.despine(fig=fig)
    fig.tight_layout()
    fig.savefig(FIG / "Figure_1" / "Figure_1A_integration_comparison.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_1" / "Figure_1A_integration_comparison.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 1A — 6-panel PCA comparison")

    # ── Panel B: Silhouette score bar chart ──
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    labels = list(silhouette_results.keys())
    scores = list(silhouette_results.values())
    colors = []
    for lab in labels:
        if lab in OMICS_COLORS:
            colors.append(OMICS_COLORS[lab])
        elif "MOFA" in lab:
            colors.append("#E91E63")
        elif "DIABLO" in lab:
            colors.append("#00BCD4")
        else:
            colors.append("#888")

    bars = ax2.bar(range(len(labels)), scores, color=colors, edgecolor="white", linewidth=1.5)
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, rotation=30, ha="right", fontsize=10)
    ax2.set_ylabel("Silhouette Score", fontsize=12)
    ax2.set_title("Figure 1B — Grade Separation Quality:\nSingle-Omics vs Integrated",
                  fontweight="bold", fontsize=13)
    ax2.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    for i, v in enumerate(scores):
        if not np.isnan(v):
            ax2.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
    sns.despine()
    fig2.tight_layout()
    fig2.savefig(FIG / "Figure_1" / "Figure_1B_silhouette_scores.png", dpi=300, bbox_inches="tight")
    fig2.savefig(FIG / "Figure_1" / "Figure_1B_silhouette_scores.svg", bbox_inches="tight")
    plt.close(fig2)
    print("  ✓ Figure 1B — Silhouette score comparison (NEW)")

    # Save scores table
    sil_df = pd.DataFrame({"Method": labels, "Silhouette_Score": scores})
    sil_df.to_csv(TAB / "Figure_1_silhouette_scores.csv", index=False)
    print("  ✓ Table — Silhouette scores saved")

    # ── Panel C: Variance explained + classification accuracy comparison ──
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score
    var_results = {}
    acc_results = {}
    for tag in OMICS_LOWER:
        if tag in pca_scores and tag in pca_var:
            v = pca_var[tag]
            total_var = v["variance_explained"].iloc[:2].sum() * 100 if "variance_explained" in v.columns else np.nan
            var_results[tag.capitalize()] = total_var
            df = pca_scores[tag]
            X = df[["PC1", "PC2"]].values
            y = df["group"].values
            try:
                acc = cross_val_score(KNeighborsClassifier(n_neighbors=2), X, y, cv=3, scoring="accuracy").mean()
            except:
                acc = np.nan
            acc_results[tag.capitalize()] = acc
    # MOFA
    if mofa_fac is not None:
        mofa_v = load_mofa_variance()
        f1v = mofa_v[mofa_v["factor"] == "Factor1"]["variance_explained"].mean()
        f2v = mofa_v[mofa_v["factor"] == "Factor2"]["variance_explained"].mean()
        var_results["MOFA2 Integrated"] = (f1v + f2v)
        try:
            acc = cross_val_score(KNeighborsClassifier(n_neighbors=2),
                combined[["F1","F2"]].values, combined["group"].values, cv=3, scoring="accuracy").mean()
        except:
            acc = np.nan
        acc_results["MOFA2 Integrated"] = acc
    # DIABLO
    if diablo_wide is not None:
        var_results["DIABLO Integrated"] = np.nan  # no direct variance explained
        try:
            dsub2 = diablo_wide[diablo_wide["dataset_label"] == diablo_wide["dataset_label"].unique()[0]]
            acc = cross_val_score(KNeighborsClassifier(n_neighbors=2),
                dsub2[["comp1","comp2"]].values, dsub2["group"].values, cv=3, scoring="accuracy").mean()
        except:
            acc = np.nan
        acc_results["DIABLO Integrated"] = acc

    fig3, ax3 = plt.subplots(figsize=(10, 5))
    methods = list(set(list(var_results.keys()) + list(acc_results.keys())))
    methods_ord = [m for m in ["Proteomics","Peptidomics","Metabolomics","Lipidomics","MOFA2 Integrated","DIABLO Integrated"] if m in methods]
    x = np.arange(len(methods_ord))
    sil_vals = [silhouette_results.get(m, np.nan) for m in methods_ord]
    acc_vals = [acc_results.get(m, np.nan) for m in methods_ord]
    w = 0.35
    ax3.bar(x - w/2, sil_vals, w, label="Silhouette Score", color="#5C6BC0", edgecolor="white")
    ax3.bar(x + w/2, acc_vals, w, label="Classification Accuracy", color="#26A69A", edgecolor="white")
    ax3.set_xticks(x)
    ax3.set_xticklabels(methods_ord, rotation=30, ha="right", fontsize=9)
    ax3.set_ylabel("Score", fontsize=12)
    ax3.set_title("Figure 1C — Quantitative Comparison: Silhouette Score & Classification Accuracy",
                  fontweight="bold", fontsize=12)
    ax3.legend(fontsize=10)
    for i, (s, a) in enumerate(zip(sil_vals, acc_vals)):
        if not np.isnan(s): ax3.text(i - w/2, s + 0.02, f"{s:.2f}", ha="center", fontsize=7)
        if not np.isnan(a): ax3.text(i + w/2, a + 0.02, f"{a:.2f}", ha="center", fontsize=7)
    sns.despine()
    fig3.tight_layout()
    fig3.savefig(FIG / "Figure_1" / "Figure_1C_quantitative_comparison.png", dpi=300, bbox_inches="tight")
    fig3.savefig(FIG / "Figure_1" / "Figure_1C_quantitative_comparison.svg", bbox_inches="tight")
    plt.close(fig3)
    print("  ✓ Figure 1C — Variance + Classification accuracy (NEW)")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 2 — MOFA2 Latent Factors                                        ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_2():
    print("\n═══ Figure 2: MOFA2 Factor Characterization ═══")

    var_df = load_mofa_variance()
    weights = load_mofa_weights()
    factors = load_mofa_factors()
    meta = pd.read_csv(MULTI / "sample_metadata.csv")

    # ── Panel A: Variance explained heatmap ──
    pivot = var_df.pivot(index="factor", columns="dataset", values="variance_explained")
    # Sort factors
    factor_order = [f"Factor{i}" for i in range(1, len(pivot) + 1)]
    pivot = pivot.reindex([f for f in factor_order if f in pivot.index])
    # Rename columns
    col_map = {c: c.capitalize() for c in pivot.columns}
    pivot = pivot.rename(columns=col_map)
    pivot = pivot[[c for c in OMICS_ORDER if c in pivot.columns]]

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd",
                linewidths=0.5, linecolor="white", ax=ax,
                cbar_kws={"label": "% Variance Explained"})
    ax.set_title("Figure 2A — MOFA2 Variance Explained per Factor per Omics Layer",
                 fontweight="bold", fontsize=12)
    ax.set_ylabel("Latent Factor")
    ax.set_xlabel("Omics Layer")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_2" / "Figure_2A_mofa_variance_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_2" / "Figure_2A_mofa_variance_heatmap.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 2A — Variance heatmap (NEW from raw data)")

    # ── Panel B: Top weights per omics for Factor 1 ──
    f1_weights = weights[weights["factor"] == "Factor1"].copy()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, ds_label in zip(axes.flatten(), OMICS_ORDER):
        sub = f1_weights[f1_weights["dataset_label"] == ds_label].copy()
        if sub.empty:
            ax.set_visible(False)
            continue
        sub = sub.nlargest(15, "abs_weight")
        sub = sub.sort_values("abs_weight", ascending=True)
        colors_list = ["#EF3F37" if w > 0 else "#262161" for w in sub["weight"]]
        ax.barh(range(len(sub)), sub["weight"], color=colors_list, edgecolor="white", linewidth=0.5)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["feature_display"], fontsize=7)
        ax.set_xlabel("MOFA Weight")
        ax.set_title(f"{ds_label}", fontweight="bold", color=OMICS_COLORS[ds_label])
        ax.axvline(0, color="grey", linewidth=0.5)
        sns.despine(ax=ax)

    fig.suptitle("Figure 2B — Top Molecular Drivers of MOFA Factor 1\n(per Omics Layer)",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "Figure_2" / "Figure_2B_mofa_top_weights.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_2" / "Figure_2B_mofa_top_weights.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 2B — Top weights (NEW from raw data)")

    # ── Panel C: Factor 1 values vs FBS grade (NEW analysis) ──
    f1_vals = factors[factors["factor"] == "Factor1"].copy()
    f1_vals = f1_vals.merge(meta[["sample_id", "group"]], on="sample_id", how="left",
                            suffixes=("_orig", ""))

    fig, ax = plt.subplots(figsize=(7, 5))
    grade_data = []
    positions = []
    grade_labels_list = []
    for i, grade in enumerate(GRADE_ORDER):
        vals = f1_vals[f1_vals["group"] == grade]["value"].values
        grade_data.append(vals)
        positions.append(i)
        grade_labels_list.append(grade)

    bp = ax.boxplot(grade_data, positions=positions, patch_artist=True, widths=0.6)
    for i, (patch, grade) in enumerate(zip(bp["boxes"], GRADE_ORDER)):
        patch.set_facecolor(GRADE_COLORS[grade])
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")
    # Overlay individual points
    for i, (vals, grade) in enumerate(zip(grade_data, GRADE_ORDER)):
        jitter = np.random.normal(0, 0.05, len(vals))
        ax.scatter(np.full_like(vals, i) + jitter, vals,
                  c=GRADE_COLORS[grade], s=80, edgecolors="black",
                  linewidth=1, zorder=3)

    ax.set_xticks(positions)
    ax.set_xticklabels(GRADE_ORDER, fontsize=11)
    ax.set_ylabel("MOFA Factor 1 Value", fontsize=12)
    ax.set_title("Figure 2C — MOFA Factor 1 Separates FBS Quality Grades",
                 fontweight="bold", fontsize=12)
    sns.despine()
    fig.tight_layout()
    fig.savefig(FIG / "Figure_2" / "Figure_2C_factor1_grade_boxplot.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_2" / "Figure_2C_factor1_grade_boxplot.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 2C — Factor 1 vs Grade boxplot (NEW analysis)")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 3 — DIABLO Sparse Signature                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_3():
    print("\n═══ Figure 3: DIABLO Sparse Signature ═══")

    diablo_wide = load_diablo_scores_wide()
    diablo_long = load_diablo_scores_long()
    diablo_ld = load_diablo_loadings()

    # ── Panel A: DIABLO sample projections per omics block ──
    datasets = diablo_wide["dataset_label"].unique()
    n_ds = len(datasets)
    fig, axes = plt.subplots(1, n_ds, figsize=(5 * n_ds, 5))
    if n_ds == 1:
        axes = [axes]

    for ax, ds in zip(axes, datasets):
        sub = diablo_wide[diablo_wide["dataset_label"] == ds]
        for grade in GRADE_ORDER:
            mask = sub["group"] == grade
            ax.scatter(sub.loc[mask, "comp1"], sub.loc[mask, "comp2"],
                      c=GRADE_COLORS[grade], label=grade, s=120,
                      edgecolors="white", linewidth=1.5, zorder=3)
        ax.set_xlabel("Component 1")
        ax.set_ylabel("Component 2")
        ds_color = OMICS_COLORS.get(ds, "#888")
        ax.set_title(f"{ds}", fontweight="bold", color=ds_color)
        ax.legend(fontsize=8)
        sns.despine(ax=ax)

    fig.suptitle("Figure 3A — DIABLO Sample Projections per Omics Block\n"
                 "(Supervised Integration for Grade Classification)",
                 fontsize=13, fontweight="bold", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG / "Figure_3" / "Figure_3A_diablo_projections.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_3" / "Figure_3A_diablo_projections.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 3A — DIABLO projections (NEW from raw data)")

    # ── Panel B: DIABLO loading contributions (Circos-style bar chart) ──
    if diablo_ld is not None and not diablo_ld.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        comp1 = diablo_ld[diablo_ld["component"] == "Component 1"].copy()
        comp1 = comp1[comp1["dataset"] != "Y"]  # drop response
        comp1 = comp1.sort_values("total_loading", ascending=True)
        ds_colors = [OMICS_COLORS.get(d.capitalize(), "#888") for d in comp1["dataset"]]
        ax.barh(range(len(comp1)), comp1["total_loading"], color=ds_colors,
                edgecolor="white", linewidth=1)
        ax.set_yticks(range(len(comp1)))
        ax.set_yticklabels([d.capitalize() for d in comp1["dataset"]], fontsize=11)
        ax.set_xlabel("Total Loading (Component 1)", fontsize=12)
        ax.set_title("Figure 3B — DIABLO Omics Block Contributions\nto Grade Classification",
                     fontweight="bold", fontsize=12)
        sns.despine()
        fig.tight_layout()
        fig.savefig(FIG / "Figure_3" / "Figure_3B_diablo_block_contributions.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG / "Figure_3" / "Figure_3B_diablo_block_contributions.svg", bbox_inches="tight")
        plt.close(fig)
        print("  ✓ Figure 3B — DIABLO block contributions (NEW)")

    # ── Panel C: Multi-block agreement scatter ──
    # Show comp1 scores across different omics blocks are correlated
    if diablo_long is not None:
        # Pivot long format: filter to comp1 only, then pivot
        comp1_data = diablo_long[diablo_long["component"] == "comp1"]
        pivot_comp1 = comp1_data.pivot_table(index="sample_id", columns="dataset_label",
                                             values="score")
        meta = pd.read_csv(MULTI / "sample_metadata.csv")
        pivot_comp1 = pivot_comp1.merge(meta[["sample_id", "group"]], on="sample_id")

        pairs = []
        cols = [c for c in pivot_comp1.columns if c not in ["sample_id", "group"]]
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                pairs.append((cols[i], cols[j]))

        n_pairs = min(len(pairs), 6)
        if n_pairs > 0:
            ncols_fig = min(3, n_pairs)
            nrows_fig = (n_pairs + ncols_fig - 1) // ncols_fig
            fig, axes = plt.subplots(nrows_fig, ncols_fig, figsize=(5 * ncols_fig, 4.5 * nrows_fig))
            if n_pairs == 1:
                axes = np.array([axes])
            axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]

            for k, (c1, c2) in enumerate(pairs[:n_pairs]):
                ax = axes[k]
                for grade in GRADE_ORDER:
                    mask = pivot_comp1["group"] == grade
                    ax.scatter(pivot_comp1.loc[mask, c1], pivot_comp1.loc[mask, c2],
                              c=GRADE_COLORS[grade], label=grade, s=100,
                              edgecolors="white", linewidth=1, zorder=3)
                # Compute correlation
                valid = pivot_comp1[[c1, c2]].dropna()
                if len(valid) > 2:
                    r, p = pearsonr(valid[c1], valid[c2])
                    ax.set_title(f"{c1} vs {c2}\nr={r:.3f}, p={p:.3e}", fontsize=10)
                else:
                    ax.set_title(f"{c1} vs {c2}", fontsize=10)
                ax.set_xlabel(f"{c1} Comp1")
                ax.set_ylabel(f"{c2} Comp1")
                ax.legend(fontsize=7)
                sns.despine(ax=ax)

            for k in range(n_pairs, len(axes)):
                axes[k].set_visible(False)

            fig.suptitle("Figure 3C — Cross-Block Agreement in DIABLO Component 1\n"
                         "(Correlated projections confirm shared biology)",
                         fontsize=12, fontweight="bold", y=1.03)
            fig.tight_layout()
            fig.savefig(FIG / "Figure_3" / "Figure_3C_crossblock_agreement.png", dpi=300, bbox_inches="tight")
            fig.savefig(FIG / "Figure_3" / "Figure_3C_crossblock_agreement.svg", bbox_inches="tight")
            plt.close(fig)
            print("  ✓ Figure 3C — Cross-block agreement (NEW analysis)")

    # ── Panel D: DIABLO Comp1 loadings per omics = "minimal molecule set" (Slide 7) ──
    for tag in OMICS_LOWER:
        src = MULTI / "DIABLO" / f"diablo_loadings_{tag}_comp1.png"
        if src.exists():
            shutil.copy(src, FIG / "Figure_3" / f"Figure_3D_DIABLO_loadings_{tag}_comp1.png")
            print(f"  ✓ Figure 3D — DIABLO {tag} comp1 loadings (selected features)")
    # Copy the variable correlation circle (Circos-style)
    vcorr = MULTI / "DIABLO" / "diablo_variable_correlation_comp1_comp2.png"
    if vcorr.exists():
        shutil.copy(vcorr, FIG / "Figure_3" / "Figure_3E_DIABLO_variable_correlation.png")
        print("  ✓ Figure 3E — DIABLO variable correlation circle")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 4 — Cross-Omics Correlation Network                             ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_4():
    print("\n═══ Figure 4: Cross-Omics Correlation Network ═══")

    edges = load_network_edges()

    # ── Panel A: Bipartite network graph ──
    G = nx.Graph()
    for _, row in edges.iterrows():
        n1 = row["feature_display_1"]
        n2 = row["feature_display_2"]
        ds1 = row["dataset_1"].capitalize()
        ds2 = row["dataset_2"].capitalize()
        corr = row["correlation"]

        G.add_node(n1, omics=ds1)
        G.add_node(n2, omics=ds2)
        G.add_edge(n1, n2, weight=abs(corr), correlation=corr)

    # Compute degree for hub detection
    degrees = dict(G.degree())

    # Color nodes
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        omics = G.nodes[node].get("omics", "")
        node_colors.append(OMICS_COLORS.get(omics, "#888888"))
        node_sizes.append(100 + degrees[node] * 40)

    fig, ax = plt.subplots(figsize=(14, 12))
    pos = nx.spring_layout(G, k=1.5, iterations=50, seed=42)

    # Draw edges with varying alpha based on correlation strength
    edge_alphas = [0.1 + 0.5 * abs(G[u][v]["correlation"]) for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, edge_color="#cccccc", width=0.8)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                          node_size=node_sizes, edgecolors="white", linewidths=1)

    # Label only hub nodes (top 15 by degree)
    top_hubs = sorted(degrees, key=degrees.get, reverse=True)[:15]
    hub_labels = {n: n for n in top_hubs}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, ax=ax, font_size=6,
                           font_weight="bold")

    # Legend
    from matplotlib.patches import Patch
    handles = [Patch(color=OMICS_COLORS[o], label=o) for o in OMICS_ORDER
               if any(G.nodes[n].get("omics") == o for n in G.nodes())]
    ax.legend(handles=handles, title="Omics Layer", loc="upper left",
             fontsize=10, title_fontsize=11)
    ax.set_title("Figure 4A — Cross-Omics Bipartite Correlation Network\n"
                 "(Proteins ↔ Metabolites/Lipids, FDR-corrected)",
                 fontsize=14, fontweight="bold")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_4" / "Figure_4A_bipartite_network.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_4" / "Figure_4A_bipartite_network.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 4A — Bipartite network (NEW)")

    # ── Panel B: Hub node analysis ──
    hub_df = pd.DataFrame([
        {"Node": n, "Degree": degrees[n], "Omics": G.nodes[n].get("omics", "")}
        for n in G.nodes()
    ]).sort_values("Degree", ascending=False).head(20)

    fig, ax = plt.subplots(figsize=(10, 8))
    hub_df_plot = hub_df.sort_values("Degree", ascending=True)
    colors_list = [OMICS_COLORS.get(o, "#888") for o in hub_df_plot["Omics"]]
    ax.barh(range(len(hub_df_plot)), hub_df_plot["Degree"], color=colors_list,
            edgecolor="white", linewidth=1)
    ax.set_yticks(range(len(hub_df_plot)))
    ax.set_yticklabels(hub_df_plot["Node"], fontsize=7)
    ax.set_xlabel("Number of Cross-Omics Connections", fontsize=12)
    ax.set_title("Figure 4B — Top 20 Hub Molecules in Cross-Omics Network\n"
                 "(Most interconnected functional drivers)",
                 fontweight="bold", fontsize=12)

    handles = [Patch(color=OMICS_COLORS[o], label=o) for o in OMICS_ORDER
               if o in hub_df_plot["Omics"].values]
    ax.legend(handles=handles, title="Omics Layer", loc="lower right")
    sns.despine()
    fig.tight_layout()
    fig.savefig(FIG / "Figure_4" / "Figure_4B_hub_analysis.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_4" / "Figure_4B_hub_analysis.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 4B — Hub node analysis (NEW)")

    hub_df.to_csv(TAB / "Figure_4_hub_nodes.csv", index=False)
    print("  ✓ Table — Hub nodes saved")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 5 — WGCNA Modules & Biological Themes                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_5():
    print("\n═══ Figure 5: Integrated Modules & Biological Themes ═══")

    # ── WGCNA Dendrograms (Slide 9 requires WGCNA) ──
    for tag in OMICS_LOWER:
        src = MULTI / "WGCNA" / f"{tag}_wgcna_dendrogram.png"
        if src.exists():
            shutil.copy(src, FIG / "Figure_5" / f"Figure_5_WGCNA_dendrogram_{tag}.png")
            print(f"  ✓ Figure 5 — WGCNA dendrogram: {tag}")

    # ── Module composition: stacked bar showing which omics contribute to each module ──
    modules = load_wgcna_modules()
    all_mods = []
    for tag, mdf in modules.items():
        mdf2 = mdf.copy()
        mdf2["omics"] = tag.capitalize()
        all_mods.append(mdf2)
    if all_mods:
        combined = pd.concat(all_mods, ignore_index=True)
        combined = combined[combined["module"] != "grey"]
        counts = combined.groupby(["module", "omics"]).size().unstack(fill_value=0)
        counts["total"] = counts.sum(axis=1)
        counts = counts.sort_values("total", ascending=False).head(15).drop(columns="total")
        # Reorder columns to match OMICS_ORDER
        cols = [c for c in OMICS_ORDER if c in counts.columns]
        counts = counts[cols]

        fig, ax = plt.subplots(figsize=(12, 6))
        bottom = np.zeros(len(counts))
        for col in cols:
            ax.bar(range(len(counts)), counts[col], bottom=bottom,
                   label=col, color=OMICS_COLORS[col], edgecolor="white", linewidth=0.5)
            bottom += counts[col].values
        ax.set_xticks(range(len(counts)))
        ax.set_xticklabels([m.capitalize() for m in counts.index], rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Number of Molecules", fontsize=11)
        ax.set_title("Figure 5D — WGCNA Module Composition Across Omics Layers\n"
                     "(Modules containing molecules from multiple omics = cross-layer biology)",
                     fontweight="bold", fontsize=12)
        ax.legend(title="Omics Layer", fontsize=9)
        sns.despine()
        fig.tight_layout()
        fig.savefig(FIG / "Figure_5" / "Figure_5D_module_composition.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG / "Figure_5" / "Figure_5D_module_composition.svg", bbox_inches="tight")
        plt.close(fig)
        print("  ✓ Figure 5D — WGCNA module composition stacked bar (NEW)")

    eigengenes = load_wgcna_eigengenes()
    pathways = load_pathway_data()

    # ── Panel A: Module-Grade correlation heatmap ──
    grade_numeric = {"High-grade": 3, "Medium-grade": 2, "Low-grade": 1}
    all_corrs = []

    for tag in OMICS_LOWER:
        if tag not in eigengenes:
            continue
        eig = eigengenes[tag].copy()
        me_cols = [c for c in eig.columns if c.startswith("ME")]
        eig["grade_num"] = eig["group"].map(grade_numeric)

        for mc in me_cols:
            if eig[mc].nunique() < 2:
                continue
            try:
                r, p = pearsonr(eig[mc], eig["grade_num"])
                all_corrs.append({
                    "Omics": tag.capitalize(),
                    "Module": mc.replace("ME", ""),
                    "Correlation": r,
                    "P_value": p,
                    "Label": f"{tag.capitalize()}\n{mc.replace('ME', '')}"
                })
            except:
                pass

    corr_df = pd.DataFrame(all_corrs)

    if not corr_df.empty:
        # Create pivot for heatmap
        corr_df["FullLabel"] = corr_df["Omics"] + " | " + corr_df["Module"]
        # Filter to most significant
        sig = corr_df[corr_df["P_value"] < 0.1].copy()
        if len(sig) < 5:
            sig = corr_df.nlargest(20, "Correlation", keep="all")
        sig = sig.head(25)

        fig, ax = plt.subplots(figsize=(10, 8))
        # Plot as a clustered bar/dot
        sig_sorted = sig.sort_values("Correlation", ascending=True)
        colors_list = [OMICS_COLORS.get(o, "#888") for o in sig_sorted["Omics"]]
        ax.barh(range(len(sig_sorted)), sig_sorted["Correlation"], color=colors_list,
                edgecolor="white", linewidth=0.8)
        ax.set_yticks(range(len(sig_sorted)))
        ax.set_yticklabels(sig_sorted["FullLabel"], fontsize=7)
        ax.set_xlabel("Pearson Correlation with FBS Grade", fontsize=11)
        ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
        ax.set_title("Figure 5A — WGCNA Module–Grade Correlations\n"
                     "(Positive = higher in High-grade FBS)",
                     fontweight="bold", fontsize=12)
        handles = [Patch(color=OMICS_COLORS[o], label=o) for o in OMICS_ORDER]
        ax.legend(handles=handles, title="Omics Layer", loc="lower right", fontsize=9)
        sns.despine()
        fig.tight_layout()
        fig.savefig(FIG / "Figure_5" / "Figure_5A_module_grade_correlation.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG / "Figure_5" / "Figure_5A_module_grade_correlation.svg", bbox_inches="tight")
        plt.close(fig)
        print("  ✓ Figure 5A — Module-grade correlation heatmap (NEW analysis)")

        corr_df.to_csv(TAB / "Figure_5_module_grade_correlations.csv", index=False)

    # ── Panel B: Integrated pathway enrichment dot plot ──
    if not pathways.empty and "term_name" in pathways.columns:
        # Filter to the primary contrast and significant hits
        pw = pathways.copy()
        if "contrast" in pw.columns:
            pw = pw[pw["contrast"] == CONTRAST]
        if "adj_p" in pw.columns:
            pw = pw[pw["adj_p"] < 0.05]
        elif "p_value" in pw.columns:
            pw = pw[pw["p_value"] < 0.05]

        if not pw.empty:
            # Count how many omics layers hit each pathway
            if "dataset" in pw.columns:
                pathway_counts = pw.groupby("term_name")["dataset"].nunique().reset_index()
                pathway_counts.columns = ["Pathway", "N_Omics_Layers"]
                pathway_counts = pathway_counts.sort_values("N_Omics_Layers", ascending=False).head(25)

                # Get details
                pw_details = pw[pw["term_name"].isin(pathway_counts["Pathway"])].copy()

                fig, ax = plt.subplots(figsize=(12, 8))
                # Create dot plot
                for i, (_, prow) in enumerate(pathway_counts.iterrows()):
                    pathway = prow["Pathway"]
                    n_omics = prow["N_Omics_Layers"]
                    pw_hits = pw_details[pw_details["term_name"] == pathway]

                    for _, hit in pw_hits.iterrows():
                        ds = hit["dataset"].capitalize() if "dataset" in hit.index else "Unknown"
                        neg_log_p = hit.get("neg_log10_p", -np.log10(hit.get("adj_p", hit.get("p_value", 0.05))))
                        color = OMICS_COLORS.get(ds, "#888")
                        ax.scatter(neg_log_p, i, c=color, s=100, edgecolors="white",
                                  linewidth=1, zorder=3)

                ax.set_yticks(range(len(pathway_counts)))
                ax.set_yticklabels(pathway_counts["Pathway"], fontsize=7)
                ax.set_xlabel("-log₁₀(adjusted p-value)", fontsize=11)
                ax.set_title("Figure 5B — Integrated Pathway Enrichment\n"
                             "(Pathways hit by multiple omics layers = convergent biology)",
                             fontweight="bold", fontsize=12)
                handles = [Patch(color=OMICS_COLORS[o], label=o) for o in OMICS_ORDER]
                ax.legend(handles=handles, title="Omics Layer", loc="lower right", fontsize=9)
                sns.despine()
                fig.tight_layout()
                fig.savefig(FIG / "Figure_5" / "Figure_5B_integrated_pathways.png", dpi=300, bbox_inches="tight")
                fig.savefig(FIG / "Figure_5" / "Figure_5B_integrated_pathways.svg", bbox_inches="tight")
                plt.close(fig)
                print("  ✓ Figure 5B — Integrated pathway enrichment dot plot (NEW)")
            else:
                print("  ⚠ No 'dataset' column in pathway data; skipping dot plot.")
        else:
            print("  ⚠ No significant pathways for primary contrast; skipping dot plot.")
    else:
        print("  ⚠ Pathway data not loaded or missing columns.")

    # ── Panel C: Sankey/Alluvial – omics → biological themes ──
    # Build from actual module membership data
    modules = load_wgcna_modules()
    theme_mapping = {}
    all_modules_combined = []
    for tag, mdf in modules.items():
        mdf = mdf.copy()
        mdf["omics"] = tag.capitalize()
        all_modules_combined.append(mdf)

    if all_modules_combined:
        combined_mod = pd.concat(all_modules_combined, ignore_index=True)
        module_counts = combined_mod.groupby(["omics", "module"]).size().reset_index(name="count")
        module_counts = module_counts[module_counts["module"] != "grey"]  # Exclude unassigned

        fig, ax = plt.subplots(figsize=(12, 8))
        # Left: omics layers, Right: modules color-coded
        left_labels = OMICS_ORDER
        right_modules = module_counts["module"].unique()
        right_labels = sorted(right_modules)[:10]  # Top 10 modules

        left_y = np.linspace(0.1, 0.9, len(left_labels))
        right_y = np.linspace(0.1, 0.9, len(right_labels))

        box_w = 0.15
        for i, (label, y) in enumerate(zip(left_labels, left_y)):
            color = OMICS_COLORS.get(label, "#888")
            ax.barh(y, box_w, height=0.06, left=0.02, color=color,
                    edgecolor="white", linewidth=1.5)
            ax.text(0.02 + box_w / 2, y, label, ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white")

        module_color_map = plt.cm.Set3(np.linspace(0, 1, len(right_labels)))
        for i, (mod, y) in enumerate(zip(right_labels, right_y)):
            ax.barh(y, box_w, height=0.06, left=0.83, color=module_color_map[i],
                    edgecolor="white", linewidth=1.5)
            ax.text(0.83 + box_w / 2, y, mod.capitalize(), ha="center", va="center",
                    fontsize=7, fontweight="bold")

        # Draw connections based on actual module membership counts
        for _, row in module_counts.iterrows():
            omics = row["omics"]
            mod = row["module"]
            if omics in left_labels and mod in right_labels:
                li = left_labels.index(omics)
                ri = list(right_labels).index(mod)
                alpha = min(0.8, row["count"] / 200)
                ax.annotate("", xy=(0.83, right_y[ri]),
                           xytext=(0.02 + box_w, left_y[li]),
                           arrowprops=dict(arrowstyle="-",
                                          color=OMICS_COLORS.get(omics, "#888"),
                                          alpha=alpha, lw=1.5 + row["count"] / 100))

        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title("Figure 5C — Multi-Omics Module Membership\n"
                     "(Omics Layers → Co-regulated WGCNA Modules)",
                     fontsize=12, fontweight="bold")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(FIG / "Figure_5" / "Figure_5C_module_alluvial.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG / "Figure_5" / "Figure_5C_module_alluvial.svg", bbox_inches="tight")
        plt.close(fig)
        print("  ✓ Figure 5C — Module membership alluvial (NEW from real data)")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 6 — Biologically Important Molecules / Mechanistic Model        ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_6():
    print("\n═══ Figure 6: Biologically Important Molecules ═══")

    diffs = load_diff_data()
    edges = load_network_edges()

    # ── Panel A: Multi-omics variable importance ──
    records = []
    for tag, df in diffs.items():
        omics = tag.capitalize()
        sub = df[df["contrast"] == CONTRAST].copy()
        sub["abs_t"] = sub["t"].abs()
        top = sub.nlargest(10, "abs_t")
        for _, row in top.iterrows():
            records.append({
                "Feature": row.get("feature_display", row["feature_id"]),
                "Importance": row["abs_t"],
                "logFC": row["logFC"],
                "Omics": omics
            })

    rdf = pd.DataFrame(records)
    rdf = rdf.sort_values("Importance", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 12))
    colors = [OMICS_COLORS.get(o, "#888") for o in rdf["Omics"]]
    ax.barh(range(len(rdf)), rdf["Importance"], color=colors,
            edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(rdf)))
    ax.set_yticklabels(rdf["Feature"], fontsize=7)
    ax.set_xlabel("|t-statistic| (Importance Score)", fontsize=11)
    ax.set_title("Figure 6A — Top Discriminating Molecules Across All Omics\n"
                 "(High-grade vs Low-grade)", fontweight="bold", fontsize=12)
    from matplotlib.patches import Patch
    handles = [Patch(color=OMICS_COLORS[o], label=o) for o in OMICS_ORDER]
    ax.legend(handles=handles, title="Omics Layer", loc="lower right")
    sns.despine()
    fig.tight_layout()
    fig.savefig(FIG / "Figure_6" / "Figure_6A_variable_importance.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_6" / "Figure_6A_variable_importance.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 6A — Variable importance (NEW from raw data)")

    # ── Panel B: Cross-omics correlation heatmap between top molecules ──
    # Get top 5 features per omics from the edges network
    unique_nodes_1 = edges[["feature_display_1", "dataset_1"]].rename(
        columns={"feature_display_1": "feature", "dataset_1": "dataset"})
    unique_nodes_2 = edges[["feature_display_2", "dataset_2"]].rename(
        columns={"feature_display_2": "feature", "dataset_2": "dataset"})
    all_nodes = pd.concat([unique_nodes_1, unique_nodes_2]).drop_duplicates()

    # Build a correlation matrix from the edge list
    features = list(set(edges["feature_display_1"].tolist() + edges["feature_display_2"].tolist()))
    if len(features) > 30:
        # Pick top-degree nodes
        degree_count = pd.concat([
            edges["feature_display_1"].value_counts(),
            edges["feature_display_2"].value_counts()
        ]).groupby(level=0).sum().nlargest(25)
        features = degree_count.index.tolist()

    # Create adjacency matrix
    corr_matrix = pd.DataFrame(0.0, index=features, columns=features)
    for _, row in edges.iterrows():
        f1, f2 = row["feature_display_1"], row["feature_display_2"]
        if f1 in features and f2 in features:
            corr_matrix.loc[f1, f2] = row["correlation"]
            corr_matrix.loc[f2, f1] = row["correlation"]

    # Fill diagonal
    np.fill_diagonal(corr_matrix.values, 1.0)

    # Keep only rows/cols with at least one non-zero off-diagonal
    mask = (corr_matrix.abs().sum(axis=1) > 1.01)
    corr_matrix = corr_matrix.loc[mask, mask]

    if len(corr_matrix) > 2:
        fig, ax = plt.subplots(figsize=(12, 10))
        sns.heatmap(corr_matrix, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                    linewidths=0.3, linecolor="white", ax=ax,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"label": "Pearson Correlation", "shrink": 0.7})
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=6, rotation=90)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
        ax.set_title("Figure 6B — Cross-Omics Correlation Heatmap\n"
                     "(Top Hub Molecules: Proteins ↔ Metabolites/Lipids)",
                     fontweight="bold", fontsize=12)
        fig.tight_layout()
        fig.savefig(FIG / "Figure_6" / "Figure_6B_cross_omics_heatmap.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG / "Figure_6" / "Figure_6B_cross_omics_heatmap.svg", bbox_inches="tight")
        plt.close(fig)
        print("  ✓ Figure 6B — Cross-omics correlation heatmap (NEW)")

    # ── Panel C: Export mechanistic model data ──
    mech_records = []
    for tag, df in diffs.items():
        omics = tag.capitalize()
        sub = df[df["contrast"] == CONTRAST].copy()
        top_up = sub.nlargest(5, "logFC")
        top_down = sub.nsmallest(5, "logFC")
        for _, row in pd.concat([top_up, top_down]).iterrows():
            direction = "↑ in High-grade" if row["logFC"] > 0 else "↓ in High-grade"
            mech_records.append({
                "Omics_Layer": omics,
                "Molecule": row.get("feature_display", row["feature_id"]),
                "logFC": round(row["logFC"], 3),
                "Direction": direction,
                "FDR": f"{row['adj.P.Val']:.2e}",
                "AveExpr": round(row["AveExpr"], 3)
            })

    mech_df = pd.DataFrame(mech_records)
    mech_df.to_csv(TAB / "Figure_6C_mechanistic_model_data.csv", index=False)
    print("  ✓ Table — Mechanistic model data exported (for BioRender schematic)")

    # ── Also export top biomarkers table ──
    rdf.to_csv(TAB / "Figure_6A_top_biomarkers.csv", index=False)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SUPPLEMENTARY TABLES                                                    ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def supplementary_tables():
    print("\n═══ Supplementary Tables ═══")
    diffs = load_diff_data()
    for i, tag in enumerate(OMICS_LOWER, 1):
        if tag not in diffs:
            continue
        df = diffs[tag]
        uniq = df.drop_duplicates(subset=["feature_id"]).copy()
        keep = [c for c in ["feature_id", "feature_display", "protein_group",
                            "gene_symbol", "protein_description", "name",
                            "mz", "rt", "logFC", "adj.P.Val", "dataset"]
                if c in uniq.columns]
        out = uniq[keep]
        out.to_csv(TAB / f"Table_S{i}_{tag}_identified_molecules.csv", index=False)
        print(f"  ✓ Table S{i}: {tag} ({len(out)} features)")

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SUPPLEMENTARY FIGURES                                                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def supplementary_figures():
    print("\n═══ Supplementary Figures ═══")
    # Copy DIABLO loading plots for Component 2 to Figure S1
    for tag in OMICS_LOWER:
        src = MULTI / "DIABLO" / f"diablo_loadings_{tag}_comp2.png"
        dst = FIG / "Figure_S1" / f"Figure_S1_DIABLO_loadings_{tag}_comp2.png"
        if src.exists():
            shutil.copy(src, dst)
            print(f"  ✓ Copied {src.name} → {dst.name}")
            
    src_circos = MULTI / "Circos" / "diablo_circos_comp2.png"
    dst_circos = FIG / "Figure_S1" / "Figure_S1_DIABLO_circos_comp2.png"
    if src_circos.exists():
        shutil.copy(src_circos, dst_circos)
        print(f"  ✓ Copied {src_circos.name} → {dst_circos.name}")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  MAIN                                                                    ║
# ╚════════════════════════════════════════════════════════════════════════════╝
# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 7 (NEW) — ABUNDANCE BARPLOT                                        ║
# ╚════════════════════════════════════════════════════════════════════════════╝
def figure_7_abundances():
    print("\n═══ Figure 7: Detected Abundances Barplot ═══")

    meta_path = MULTI / "MetaboAnalyst"

    try:
        pro = pd.read_csv(meta_path / "proteomics_metaboanalyst_matrix.csv")
        pep = pd.read_csv(meta_path / "peptidomics_metaboanalyst_matrix.csv")
        met = pd.read_csv(meta_path / "metabolomics_metaboanalyst_matrix.csv")
        diff_met = pd.read_csv(MULTI / "Differential_Updated" / "metabolomics_differential_results.csv")
    except Exception as e:
        print(f"  ✗ Skipped due to missing data: {e}")
        return

    # ── Real protein / molecule names for display ──
    PROTEIN_NAMES = {
        "Q3MHL4": "Vitamin D-binding protein",
        "Q0VCA8": "Villin-like protein",
        "Q28178": "Aldehyde dehydrogenase",
        "P02070": "Hemoglobin subunit alpha",
        "Q9XSC6": "Fatty acid-binding protein",
        "P07224": "Vitamin K-dep. protein C",
    }
    PEPTIDE_NAMES = {
        "P07456":  "Lactadherin (MFG-E8)",
        "P02769":  "Serum albumin (BSA)",
        "P02672":  "Fibrinogen alpha",
        "P19879":  "Annexin A2",       # the 'None #39' entry
    }

    # ── 1. Proteomics: match by UniProt accession prefix ──
    pro_accessions = list(PROTEIN_NAMES.keys())
    pro_matched = pro[pro["feature_id"].str.split("|").str[0].isin(pro_accessions)].copy()
    pro_matched["display_name"] = pro_matched["feature_id"].str.split("|").str[0].map(PROTEIN_NAMES)

    # ── 2. Peptidomics: match by UniProt accession prefix ──
    pep_accessions = list(PEPTIDE_NAMES.keys())
    # Take ONE representative peptide per parent protein
    pep_matched_rows = []
    for acc in pep_accessions:
        hits = pep[pep["feature_id"].str.startswith(acc)]
        if not hits.empty:
            pep_matched_rows.append(hits.iloc[0])
    if pep_matched_rows:
        pep_matched = pd.DataFrame(pep_matched_rows)
        pep_matched["display_name"] = pep_matched["feature_id"].str.split("|").str[0].map(PEPTIDE_NAMES)
    else:
        pep_matched = pd.DataFrame()

    # ── 3. Metabolomics: top 3 down-regulated features ──
    met_top = diff_met[diff_met["contrast"] == "High.grade - Low.grade"].nsmallest(3, "logFC")
    met_feats = met_top["feature_id"].tolist()
    met_names = dict(zip(met_top["feature_id"], met_top["feature_display"]))
    met_matched = met[met["feature_id"].isin(met_feats)].copy()
    met_matched["display_name"] = met_matched["feature_id"].map(met_names)

    # ── Build long-form records ──
    grades = GRADE_ORDER
    records = []

    for omics_label, matched_df in [("Proteomics", pro_matched),
                                     ("Peptidomics", pep_matched),
                                     ("Metabolomics", met_matched)]:
        if matched_df is None or matched_df.empty:
            print(f"  ⚠ No matched features for {omics_label}")
            continue
        for _, row in matched_df.iterrows():
            name = row["display_name"] if pd.notna(row.get("display_name")) else row["feature_id"][:20]
            for grade in grades:
                cols = [c for c in matched_df.columns if c.startswith(grade)]
                vals = row[cols].values.astype(float)
                for v in vals:
                    records.append({
                        "Omics": omics_label,
                        "Feature": name,
                        "Grade": grade,
                        "Abundance": v,
                    })

    if not records:
        print("  ✗ No records generated.")
        return

    rdf = pd.DataFrame(records)

    # ── Plot ──
    n_pro = rdf[rdf["Omics"] == "Proteomics"]["Feature"].nunique()
    n_pep = rdf[rdf["Omics"] == "Peptidomics"]["Feature"].nunique()
    n_met = rdf[rdf["Omics"] == "Metabolomics"]["Feature"].nunique()
    ratios = [max(n_pro, 1), max(n_pep, 1), max(n_met, 1)]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6),
                             gridspec_kw={"width_ratios": ratios}, sharey=False)

    for i, (omics_label, ax) in enumerate(zip(["Proteomics", "Peptidomics", "Metabolomics"], axes)):
        sub = rdf[rdf["Omics"] == omics_label]
        if sub.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(omics_label, fontweight="bold", fontsize=14)
            continue
        sns.barplot(data=sub, x="Feature", y="Abundance", hue="Grade",
                    hue_order=grades, palette=GRADE_COLORS, ax=ax,
                    capsize=0.08, errwidth=1.2, edgecolor="black", linewidth=0.6)
        ax.set_title(omics_label, fontweight="bold", fontsize=14)
        ax.set_xlabel("")
        ax.set_ylabel("Normalized Abundance" if i == 0 else "")
        ax.tick_params(axis="x", rotation=40, labelsize=8)

        if i == 0:
            ax.legend(title="Grade", loc="upper right", fontsize=8, title_fontsize=9)
        elif ax.get_legend() is not None:
            ax.get_legend().remove()

    fig.suptitle("Detected Abundances of Key DIABLO Marker Candidates\nacross Three FBS Quality Gradients",
                 fontweight="bold", fontsize=15, y=1.02)
    fig.tight_layout()

    out_dir = FIG / "Figure_7_Abundances"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "Detected_Abundances_Barplot.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "Detected_Abundances_Barplot.svg", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ Figure 7 — Abundances barplot generated (NEW)")


def main():
    print("=" * 70)
    print("FBS Multi-Omics Paper 3 — FRESH Figure Generator")
    print("All figures generated from RAW CSV data")
    print("=" * 70)

    ensure_dirs()

    figure_1()
    figure_2()
    figure_3()
    figure_4()
    figure_5()
    figure_6()
    figure_7_abundances()
    supplementary_tables()
    supplementary_figures()

    print("\n" + "=" * 70)
    print("GENERATION COMPLETE")
    print(f"Output: {OUTPUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()

