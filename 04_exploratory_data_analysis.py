#!/usr/bin/env python3
"""
04_exploratory_data_analysis.py
Exploratory data analysis of Nordic pollen data (Greenland excluded)
High-quality outputs (300 dpi). Period bar plot removed.
Creates a combined 3-map TIFF with vertical subplots and enhanced visuals.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from mpl_toolkits.basemap import Basemap
import warnings
warnings.filterwarnings('ignore')

# --- Global style & larger fonts ---
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 20,
    "axes.labelsize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 22
})

def load_data():
    print("="*70)
    print("NORDIC POLLEN DATA - EXPLORATORY DATA ANALYSIS")
    print("="*70)

    output_dir = Path("04_exploratory_data_analysis")
    output_dir.mkdir(exist_ok=True)

    df = pd.read_excel("03_nordic_analysis/nordic_pollen_data.xlsx")
    print(f"\nDataset loaded: {df.shape[0]:,} records, {df.shape[1]} columns")
    return df, output_dir

def basic_statistics(df, output_dir):
    print("\n" + "="*70)
    print("1. BASIC STATISTICS & DATA QUALITY")
    print("-"*70)

    stats_dict = {
        'Total Records': len(df),
        'Unique Sites': df['site_name'].nunique() if 'site_name' in df else np.nan,
        'Unique Taxa': df['taxon'].nunique() if 'taxon' in df else np.nan,
        'Unique Samples': df['sample_id'].nunique() if 'sample_id' in df else np.nan,
        'Unique Datasets': df['dataset_id'].nunique() if 'dataset_id' in df else np.nan,
    }

    missing_stats = df.isnull().sum()
    missing_pct = (missing_stats / len(df) * 100).round(2)

    print("\nBasic Statistics:")
    for k, v in stats_dict.items():
        print(f"  {k:20s}: {int(v):,}" if pd.notna(v) else f"  {k:20s}: N/A")

    print("\nMissing Values:")
    for col in df.columns:
        if missing_stats[col] > 0:
            print(f"  {col:20s}: {missing_stats[col]:,} ({missing_pct[col]:.1f}%)")

    with pd.ExcelWriter(output_dir / 'data_quality_report.xlsx') as writer:
        pd.DataFrame(list(stats_dict.items()), columns=['Metric', 'Value']).to_excel(
            writer, sheet_name='Summary', index=False)
        pd.DataFrame({'Column': missing_stats.index, 'Missing_Count': missing_stats.values,
                      'Missing_Pct': missing_pct.values}).to_excel(
            writer, sheet_name='Missing_Values', index=False)

    return stats_dict

def temporal_analysis(df, output_dir):
    """ONLY the age histogram with beautiful colors."""
    print("\n" + "="*70)
    print("2. TEMPORAL ANALYSIS (Histogram only)")
    print("-"*70)

    # Valid ages only
    df_valid = df[df['age_BP'].ge(0)].copy()
    neg_n = len(df) - len(df_valid)
    if neg_n > 0:
        print(f"Excluded {neg_n:,} records with negative age_BP")

    # Histogram with beautiful gradient colors
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create histogram with beautiful teal/turquoise color
    n, bins, patches = ax.hist(df_valid['age_BP'], bins=60, edgecolor='white', 
                                alpha=0.85, linewidth=1.2)
    
    # Apply gradient colors to bars (fixed for newer matplotlib)
    from matplotlib import cm
    cmap = cm.viridis
    norm = plt.Normalize(vmin=bins.min(), vmax=bins.max())
    for i, (patch, bin_val) in enumerate(zip(patches, bins[:-1])):
        patch.set_facecolor(cmap(norm(bin_val)))
    
    ax.set_xlabel('Age (BP)')
    ax.set_ylabel('Number of Records')
    ax.set_title('Age Distribution of Nordic Pollen Records')
    ax.axvline(11_700, color='#FF6B6B', linestyle='--', linewidth=3, 
               label='Holocene boundary (11,700 BP)')
    
    # Legend outside the plot box
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, 
              frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'temporal_full_age_distribution.tif', dpi=300, 
                format='tiff', bbox_inches='tight')
    plt.close()

def _basemap_nordic(ax):
    """Nordic-only map (excludes Greenland)."""
    m = Basemap(projection='cyl',
                llcrnrlat=50, urcrnrlat=82,
                llcrnrlon=-30, urcrnrlon=40,
                resolution='i', ax=ax)
    m.shadedrelief(scale=1)
    m.drawcoastlines(color='dimgray', linewidth=0.6)
    m.drawcountries(color='dimgray', linewidth=0.6)
    parallels = np.arange(50, 83, 5)
    meridians = np.arange(-30, 41, 10)
    m.drawparallels(parallels, labels=[1,0,0,0], fontsize=12, linewidth=0.4, color='gray')
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=12, linewidth=0.4, color='gray')
    return m

def _compute_site_locs(df):
    site_locs = df.groupby('site_name').agg({
        'latitude': 'first',
        'longitude': 'first',
        'sample_id': 'nunique',
        'taxon': 'nunique',
        'age_BP': ['mean', 'min', 'max']
    }).reset_index()
    site_locs.columns = ['site_name', 'latitude', 'longitude', 'n_samples',
                         'n_taxa', 'mean_age', 'min_age', 'max_age']
    return site_locs

def spatial_analysis_basemap(df, output_dir):
    """Create individual maps AND a combined vertical 3-map figure."""
    print("\n" + "="*70)
    print("3. SPATIAL ANALYSIS - BASEMAP VISUALIZATION (Nordic only)")
    print("-"*70)

    site_locs = _compute_site_locs(df)

    # ---- Individual map 1: Samples ----
    fig, ax = plt.subplots(figsize=(20, 14))
    m = _basemap_nordic(ax)
    vals = site_locs['n_samples'].to_numpy()
    vmax = np.percentile(vals, 95) if len(vals) > 5 else (vals.max() if len(vals) else 1)
    norm_size = np.power(vals / max(vals.max(), 1), 0.5) if len(vals) else np.array([0])

    sc = m.scatter(site_locs['longitude'].values, site_locs['latitude'].values,
                   latlon=True, c=vals, s=200 + 300*norm_size,
                   cmap='plasma', alpha=0.9, edgecolors='white', linewidth=1.6,
                   vmin=0, vmax=vmax, zorder=5)
    cbar = plt.colorbar(sc, orientation='vertical', pad=0.02, shrink=0.75, aspect=24)
    cbar.set_label('Number of Samples per Site')
    ax.set_title('Nordic Pollen Sites — Sample Distribution')
    plt.tight_layout()
    plt.savefig(output_dir / 'spatial_basemap_samples.tif', dpi=300, format='tiff', bbox_inches='tight')
    plt.close()

    # ---- Individual map 2: Diversity ----
    fig, ax = plt.subplots(figsize=(20, 14))
    m = _basemap_nordic(ax)
    vals = site_locs['n_taxa'].to_numpy()
    vmax = np.percentile(vals, 95) if len(vals) > 5 else (vals.max() if len(vals) else 1)

    sc = m.scatter(site_locs['longitude'].values, site_locs['latitude'].values,
                   latlon=True, c=vals, s=260,
                   cmap='viridis', alpha=0.9, edgecolors='white', linewidth=1.6,
                   vmin=0, vmax=vmax, zorder=5)
    cbar = plt.colorbar(sc, orientation='vertical', pad=0.02, shrink=0.75, aspect=24)
    cbar.set_label('Number of Taxa per Site')
    ax.set_title('Nordic Pollen Sites — Taxonomic Diversity')
    plt.tight_layout()
    plt.savefig(output_dir / 'spatial_basemap_diversity.tif', dpi=300, format='tiff', bbox_inches='tight')
    plt.close()

    # ---- Individual map 3: Temporal coverage ----
    fig, ax = plt.subplots(figsize=(20, 14))
    m = _basemap_nordic(ax)
    vals = (site_locs['max_age'] - site_locs['min_age']).to_numpy()
    vmax = np.percentile(vals, 95) if len(vals) > 5 else (vals.max() if len(vals) else 1)

    sc = m.scatter(site_locs['longitude'].values, site_locs['latitude'].values,
                   latlon=True, c=vals, s=260,
                   cmap='coolwarm', alpha=0.9, edgecolors='white', linewidth=1.6,
                   vmin=0, vmax=vmax, zorder=5)
    cbar = plt.colorbar(sc, orientation='vertical', pad=0.02, shrink=0.75, aspect=24)
    cbar.set_label('Temporal Coverage (years)')
    ax.set_title('Nordic Pollen Sites — Temporal Coverage')
    plt.tight_layout()
    plt.savefig(output_dir / 'spatial_basemap_temporal.tif', dpi=300, format='tiff', bbox_inches='tight')
    plt.close()

    # ---- Combined VERTICAL 3-map figure ----
    fig, axes = plt.subplots(3, 1, figsize=(22, 42))
    
    # Map 1: Samples (top)
    ax_top = axes[0]
    m = _basemap_nordic(ax_top)
    vals = site_locs['n_samples'].to_numpy()
    vmax = np.percentile(vals, 95) if len(vals) > 5 else (vals.max() if len(vals) else 1)
    norm_size = np.power(vals / max(vals.max(), 1), 0.5) if len(vals) else np.array([0])
    sc = m.scatter(site_locs['longitude'].values, site_locs['latitude'].values,
                   latlon=True, c=vals, s=200 + 300*norm_size,
                   cmap='plasma', alpha=0.9, edgecolors='white', linewidth=1.6,
                   vmin=0, vmax=vmax, zorder=5)
    cbar = fig.colorbar(sc, ax=ax_top, orientation='vertical', pad=0.02, 
                        shrink=0.85, aspect=28)
    cbar.set_label('Number of Samples per Site', fontsize=16)
    ax_top.set_title('(A) Sample Distribution', fontsize=22, pad=20)

    # Map 2: Diversity (middle)
    ax_mid = axes[1]
    m = _basemap_nordic(ax_mid)
    vals = site_locs['n_taxa'].to_numpy()
    vmax = np.percentile(vals, 95) if len(vals) > 5 else (vals.max() if len(vals) else 1)
    sc = m.scatter(site_locs['longitude'].values, site_locs['latitude'].values,
                   latlon=True, c=vals, s=260,
                   cmap='viridis', alpha=0.9, edgecolors='white', linewidth=1.6,
                   vmin=0, vmax=vmax, zorder=5)
    cbar = fig.colorbar(sc, ax=ax_mid, orientation='vertical', pad=0.02, 
                        shrink=0.85, aspect=28)
    cbar.set_label('Number of Taxa per Site', fontsize=16)
    ax_mid.set_title('(B) Taxonomic Diversity', fontsize=22, pad=20)

    # Map 3: Temporal coverage (bottom)
    ax_bot = axes[2]
    m = _basemap_nordic(ax_bot)
    vals = (site_locs['max_age'] - site_locs['min_age']).to_numpy()
    vmax = np.percentile(vals, 95) if len(vals) > 5 else (vals.max() if len(vals) else 1)
    sc = m.scatter(site_locs['longitude'].values, site_locs['latitude'].values,
                   latlon=True, c=vals, s=260,
                   cmap='coolwarm', alpha=0.9, edgecolors='white', linewidth=1.6,
                   vmin=0, vmax=vmax, zorder=5)
    cbar = fig.colorbar(sc, ax=ax_bot, orientation='vertical', pad=0.02, 
                        shrink=0.85, aspect=28)
    cbar.set_label('Temporal Coverage (years)', fontsize=16)
    ax_bot.set_title('(C) Temporal Coverage', fontsize=22, pad=20)

    fig.suptitle("Nordic Pollen Sites — Spatial Patterns", fontsize=26, y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / 'spatial_basemap_ALL_IN_ONE_VERTICAL.tif', 
                dpi=300, format='tiff', bbox_inches='tight')
    plt.close()

    print("\nCreated Basemap visualizations (Nordic only):")
    print("  - spatial_basemap_samples.tif")
    print("  - spatial_basemap_diversity.tif")
    print("  - spatial_basemap_temporal.tif")
    print("  - spatial_basemap_ALL_IN_ONE_VERTICAL.tif (vertical layout with enhanced depth)")

    return site_locs

def site_analysis(df, output_dir):
    print("\n" + "="*70)
    print("4. SITE-LEVEL ANALYSIS")
    print("-"*70)

    site_stats = df.groupby('site_name').agg({
        'sample_id': 'nunique',
        'taxon': 'nunique',
        'age_BP': ['min', 'max', 'mean'],
        'latitude': 'first',
        'longitude': 'first',
        'count': 'sum' if 'count' in df.columns else 'size'
    }).round(2)

    site_stats.columns = ['n_samples', 'n_taxa', 'min_age', 'max_age', 'mean_age',
                          'latitude', 'longitude', 'total_count']
    site_stats['age_range'] = site_stats['max_age'] - site_stats['min_age']

    site_stats.sort_values('n_samples', ascending=False).to_excel(output_dir / 'site_statistics.xlsx')
    return site_stats

def generate_summary_report(df, output_dir, _all_stats):
    print("\n" + "="*70)
    print("5. GENERATING SUMMARY REPORT")
    print("-"*70)

    negative_age_count = (df['age_BP'] < 0).sum()
    report_path = output_dir / 'EDA_SUMMARY_REPORT.txt'

    with open(report_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("NORDIC POLLEN DATA - EXPLORATORY DATA ANALYSIS SUMMARY\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total records loaded: {len(df):,}\n")
        f.write(f"Records with negative age_BP (excluded): {negative_age_count:,}\n")
        f.write(f"Valid records analyzed: {len(df[df['age_BP'] >= 0]):,}\n")
        f.write(f"Unique sites: {df['site_name'].nunique() if 'site_name' in df else 'N/A'}\n")
        f.write(f"Unique taxa: {df['taxon'].nunique() if 'taxon' in df else 'N/A'}\n")
        if (df['age_BP'] >= 0).any():
            f.write(f"Temporal range: {df[df['age_BP'] >= 0]['age_BP'].min():.0f} to "
                    f"{df[df['age_BP'] >= 0]['age_BP'].max():.0f} BP\n")
        f.write(f"Spatial extent: {df['latitude'].min():.1f}-{df['latitude'].max():.1f}°N, "
                f"{df['longitude'].min():.1f}-{df['longitude'].max():.1f}°E\n\n")
        f.write("All figures saved at 300 dpi.\n")

    print(f"\nSummary report saved to: {report_path}")
    print("\nEDA COMPLETE! All results saved to: 04_exploratory_data_analysis/")

def main():
    df, output_dir = load_data()
    basic_stats = basic_statistics(df, output_dir)
    temporal_analysis(df, output_dir)
    site_locs = spatial_analysis_basemap(df, output_dir)
    site_stats = site_analysis(df, output_dir)
    generate_summary_report(df, output_dir, {
        'basic': basic_stats,
        'spatial': site_locs,
        'sites': site_stats
    })
    print("\n" + "="*70)
    print("EXPLORATORY DATA ANALYSIS COMPLETE")
    print("="*70)

if __name__ == "__main__":
    main()