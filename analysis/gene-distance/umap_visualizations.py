"""
UMAP visualizations for gene-distance analysis.

Loads processed RNA/protein data, builds a shared UMAP space (fit on RNA,
project protein), and creates publication-quality plots showing the
RNA-protein distance for each gene.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from umap import UMAP
import os
from pathlib import Path

# Configuration
DATA_DIR = Path(__file__).parent.parent.parent / "data" / "real"
OUTPUT_DIR = Path(__file__).parent / "figs"
OUTPUT_DIR.mkdir(exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 5)
plt.rcParams['font.size'] = 10


def load_variety_data(variety_name):
    """Load RNA and protein data for a variety (Cadenza or Norin)."""
    lower_name = variety_name.lower()

    # Load counts/intensities
    rna = pd.read_csv(f"{DATA_DIR}/{lower_name}-rnaseq.csv", index_col=0)
    prot = pd.read_csv(f"{DATA_DIR}/{lower_name}-prot.csv", index_col=0)

    # Load metadata
    meta = pd.read_csv(f"{DATA_DIR}/{lower_name}_metadata.csv", index_col=0)

    # Load gene mapping
    gmap = pd.read_csv(f"{DATA_DIR}/{lower_name}_protein_gene_mapping.csv")

    return {"rna": rna, "prot": prot, "meta": meta, "gmap": gmap}


def cell_means(data, meta, condition_col="condition"):
    """Compute condition-level means from replicate data.

    Parameters
    ----------
    data : pd.DataFrame
        Gene x sample matrix (samples as columns)
    meta : pd.DataFrame
        Sample metadata with condition column
    condition_col : str
        Column name defining condition groups

    Returns
    -------
    pd.DataFrame
        Gene x condition matrix of means
    """
    common_idx = data.columns.intersection(meta.index)
    data_common = data[common_idx]
    meta_common = meta.loc[common_idx]

    grouped = meta_common.groupby(condition_col).groups
    means = pd.DataFrame(index=data_common.index)

    for cond, samples in grouped.items():
        means[cond] = data_common[samples].mean(axis=1)

    return means


def row_center(M):
    """Subtract each row's mean (center by gene)."""
    return M - M.mean(axis=1, keepdims=True)


def fit_col_scaler(M):
    """Fit per-column (per-condition) z-score standardization on RNA only."""
    mu = M.mean(axis=0)
    sdv = M.std(axis=0)
    sdv[sdv < 1e-8] = 1
    return mu, sdv


def apply_col_scaler(M, mu, sdv):
    """Apply fitted column-wise standardization."""
    return (M - mu) / sdv


def build_shared_space_umap(R_cond, P_cond, n_neighbors=15, min_dist=0.1, seed=42):
    """Build shared UMAP space: fit on RNA, project protein.

    Parameters
    ----------
    R_cond : pd.DataFrame
        Gene x condition matrix (RNA)
    P_cond : pd.DataFrame
        Gene x condition matrix (protein, same genes/order as R_cond)
    n_neighbors : int
        UMAP n_neighbors parameter
    min_dist : float
        UMAP min_dist parameter
    seed : int
        Random seed

    Returns
    -------
    dict with keys:
        - Rs, Ps: standardized matrices
        - umap_fit: fitted UMAP model
        - Er, Ep: UMAP embeddings (RNA, protein)
        - R_cond, P_cond: original condition-mean matrices
    """
    # Row-center (compare response shape, not level)
    R_centered = row_center(R_cond)
    P_centered = row_center(P_cond)

    # Fit scaler on RNA only, apply to both
    mu, sdv = fit_col_scaler(R_centered.values)
    Rs = apply_col_scaler(R_centered.values, mu, sdv)
    Ps = apply_col_scaler(P_centered.values, mu, sdv)

    Rs_df = pd.DataFrame(Rs, index=R_cond.index, columns=R_cond.columns)
    Ps_df = pd.DataFrame(Ps, index=P_cond.index, columns=P_cond.columns)

    # Fit UMAP on RNA, transform protein
    nn = max(2, min(15, Rs.shape[0] // 4))
    umap_fit = UMAP(
        n_components=2,
        n_neighbors=nn,
        min_dist=min_dist,
        metric='euclidean',
        random_state=seed
    )
    Er = umap_fit.fit_transform(Rs)
    Ep = umap_fit.transform(Ps)

    Er_df = pd.DataFrame(Er, index=R_cond.index, columns=['UMAP1', 'UMAP2'])
    Ep_df = pd.DataFrame(Ep, index=P_cond.index, columns=['UMAP1', 'UMAP2'])

    return {
        'Rs': Rs_df, 'Ps': Ps_df,
        'Er': Er_df, 'Ep': Ep_df,
        'umap_fit': umap_fit,
        'R_cond': R_cond, 'P_cond': P_cond
    }


def compute_distances(space_dict):
    """Compute distances between RNA and protein in different spaces.

    Parameters
    ----------
    space_dict : dict
        Output from build_shared_space_umap

    Returns
    -------
    pd.DataFrame
        Genes x distance metrics
    """
    Rs = space_dict['Rs'].values
    Ps = space_dict['Ps'].values
    Er = space_dict['Er'].values
    Ep = space_dict['Ep'].values

    d_full = np.sqrt(((Rs - Ps) ** 2).sum(axis=1))
    d_umap = np.sqrt(((Er - Ep) ** 2).sum(axis=1))

    dist_df = pd.DataFrame({
        'gene_id': space_dict['Rs'].index,
        'd_full': d_full,
        'd_umap': d_umap
    })

    return dist_df


def plot_shared_umap(space_dict, dist_df, variety_name, output_path):
    """Create UMAP scatter plot with RNA and protein points.

    Shows RNA and protein positions in shared UMAP space, with points
    colored by distance between them.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    Er = space_dict['Er'].values
    Ep = space_dict['Ep'].values
    distances = dist_df['d_umap'].values

    # Left panel: colored by distance magnitude
    scatter = ax1.scatter(Er[:, 0], Er[:, 1], c=distances, cmap='RdYlBu_r',
                         s=80, alpha=0.6, label='RNA', edgecolors='black', linewidth=0.5)
    ax1.scatter(Ep[:, 0], Ep[:, 1], c=distances, cmap='RdYlBu_r',
               s=80, alpha=0.6, marker='^', label='Protein', edgecolors='black', linewidth=0.5)

    # Draw lines connecting RNA-protein pairs
    for i in range(len(Er)):
        ax1.plot([Er[i, 0], Ep[i, 0]], [Er[i, 1], Ep[i, 1]],
                'k-', alpha=0.1, linewidth=0.5)

    ax1.set_xlabel('UMAP1')
    ax1.set_ylabel('UMAP2')
    ax1.set_title(f'{variety_name}: RNA-Protein Shared UMAP Space\n(colored by distance)')
    ax1.legend()
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('RNA-Protein Distance')

    # Right panel: RNA (circles) and Protein (triangles) separately colored
    ax2.scatter(Er[:, 0], Er[:, 1], s=80, alpha=0.7, c='#2E86AB',
               label='RNA', edgecolors='black', linewidth=0.5)
    ax2.scatter(Ep[:, 0], Ep[:, 1], s=80, alpha=0.7, c='#A23B72', marker='^',
               label='Protein', edgecolors='black', linewidth=0.5)

    ax2.set_xlabel('UMAP1')
    ax2.set_ylabel('UMAP2')
    ax2.set_title(f'{variety_name}: RNA vs Protein Positions')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_distance_histogram(dist_df, variety_name, output_path):
    """Create histogram of RNA-protein distances."""
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(dist_df['d_full'], bins=30, alpha=0.7, color='#2E86AB', edgecolor='black')
    ax.axvline(dist_df['d_full'].mean(), color='red', linestyle='--', linewidth=2, label='Mean')
    ax.axvline(dist_df['d_full'].median(), color='orange', linestyle='--', linewidth=2, label='Median')

    ax.set_xlabel('RNA-Protein Distance (full space)')
    ax.set_ylabel('Count')
    ax.set_title(f'{variety_name}: Distribution of Gene-Level RNA-Protein Distances')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_condition_means(space_dict, variety_name, output_path):
    """Create heatmaps of condition-mean profiles."""
    R_cond = space_dict['R_cond']
    P_cond = space_dict['P_cond']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Normalize for visualization
    R_norm = (R_cond - R_cond.mean()) / R_cond.std()
    P_norm = (P_cond - P_cond.mean()) / P_cond.std()

    # RNA heatmap
    sns.heatmap(R_norm.iloc[:min(50, len(R_norm))], cmap='RdBu_r', center=0,
               cbar_kws={'label': 'Z-score'}, ax=ax1, yticklabels=False)
    ax1.set_title(f'{variety_name}: RNA Condition-Mean Profiles (top 50 genes)')
    ax1.set_xlabel('Condition')

    # Protein heatmap
    sns.heatmap(P_norm.iloc[:min(50, len(P_norm))], cmap='RdBu_r', center=0,
               cbar_kws={'label': 'Z-score'}, ax=ax2, yticklabels=False)
    ax2.set_title(f'{variety_name}: Protein Condition-Mean Profiles (top 50 genes)')
    ax2.set_xlabel('Condition')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_top_genes(dist_df, space_dict, variety_name, output_path, n_top=10):
    """Highlight top/bottom genes by distance."""
    fig, ax = plt.subplots(figsize=(12, 8))

    Er = space_dict['Er'].values
    Ep = space_dict['Ep'].values
    distances = dist_df['d_umap'].values

    # All genes in light gray
    ax.scatter(Er[:, 0], Er[:, 1], s=30, alpha=0.2, c='gray', label='RNA (background)')
    ax.scatter(Ep[:, 0], Ep[:, 1], s=30, alpha=0.2, c='gray', marker='^', label='Protein (background)')

    # Top N most discordant (high distance)
    top_idx = np.argsort(-distances)[:n_top]
    ax.scatter(Er[top_idx, 0], Er[top_idx, 1], s=150, alpha=0.8, c='#D1495B',
              edgecolors='black', linewidth=1.5, label=f'RNA (top {n_top} discordant)')
    ax.scatter(Ep[top_idx, 0], Ep[top_idx, 1], s=150, alpha=0.8, c='#D1495B', marker='^',
              edgecolors='black', linewidth=1.5, label=f'Protein (top {n_top} discordant)')

    # Draw connecting lines for top genes
    for i in top_idx:
        ax.plot([Er[i, 0], Ep[i, 0]], [Er[i, 1], Ep[i, 1]],
               '#D1495B', alpha=0.5, linewidth=1.5)

    # Bottom N most concordant (low distance)
    bottom_idx = np.argsort(distances)[:n_top]
    ax.scatter(Er[bottom_idx, 0], Er[bottom_idx, 1], s=150, alpha=0.8, c='#2C6FBB',
              edgecolors='black', linewidth=1.5, label=f'RNA (bottom {n_top} concordant)')
    ax.scatter(Ep[bottom_idx, 0], Ep[bottom_idx, 1], s=150, alpha=0.8, c='#2C6FBB', marker='^',
              edgecolors='black', linewidth=1.5, label=f'Protein (bottom {n_top} concordant)')

    # Draw connecting lines for bottom genes
    for i in bottom_idx:
        ax.plot([Er[i, 0], Ep[i, 0]], [Er[i, 1], Ep[i, 1]],
               '#2C6FBB', alpha=0.5, linewidth=1.5)

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title(f'{variety_name}: Top {n_top} Most/Least Discordant Genes')
    ax.legend(loc='best', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def main():
    """Main pipeline: load data, build shared spaces, create visualizations."""

    for variety_name in ['Cadenza', 'Norin']:
        print(f"\n{'='*60}")
        print(f"Processing {variety_name}")
        print('='*60)

        # Load data
        data = load_variety_data(variety_name)
        print(f"RNA shape: {data['rna'].shape}")
        print(f"Protein shape: {data['prot'].shape}")

        # Find common genes (1:1 mapping)
        common_genes = set(data['rna'].index) & set(data['prot'].index)
        print(f"Common genes: {len(common_genes)}")

        # Compute condition means
        R_cond = cell_means(data['rna'], data['meta'], condition_col='condition')
        P_cond = cell_means(data['prot'], data['meta'], condition_col='condition')

        # Keep only common genes
        R_cond = R_cond.loc[common_genes]
        P_cond = P_cond.loc[common_genes]

        print(f"Condition-mean matrices: {R_cond.shape}")
        print(f"Conditions: {list(R_cond.columns)}")

        # Build shared UMAP space
        space_dict = build_shared_space_umap(R_cond, P_cond)
        print("Built shared UMAP space")

        # Compute distances
        dist_df = compute_distances(space_dict)
        print(f"Distance stats:")
        print(f"  Full space: {dist_df['d_full'].describe()}")
        print(f"  UMAP: {dist_df['d_umap'].describe()}")

        # Create visualizations
        prefix = variety_name.lower()

        plot_shared_umap(
            space_dict, dist_df, variety_name,
            OUTPUT_DIR / f"{prefix}_shared_umap.png"
        )

        plot_distance_histogram(
            dist_df, variety_name,
            OUTPUT_DIR / f"{prefix}_distance_histogram.png"
        )

        plot_condition_means(
            space_dict, variety_name,
            OUTPUT_DIR / f"{prefix}_condition_means.png"
        )

        plot_top_genes(
            dist_df, space_dict, variety_name,
            OUTPUT_DIR / f"{prefix}_top_genes.png",
            n_top=10
        )

        # Save distance table
        dist_df_sorted = dist_df.sort_values('d_full', ascending=False)
        dist_df_sorted.to_csv(
            OUTPUT_DIR / f"{prefix}_gene_distances.csv",
            index=False
        )
        print(f"Saved distance table: {prefix}_gene_distances.csv")


if __name__ == '__main__':
    main()
