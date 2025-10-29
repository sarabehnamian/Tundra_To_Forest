## Tundra To Forest – Pollen-Based Paleovegetation Workflow

This repository contains a reproducible pipeline to download, clean, analyze, and visualize paleoecological pollen records to study tundra-to-forest transitions. It automates data retrieval from Neotoma, merges and fixes records, performs Nordic-focused analyses, generates exploratory data summaries, and computes biotic time series metrics and plots.

### Repository structure

```
D:\Tundra_To_Forest\
  00_neotoma_download_extract_geo.py
  01_merge_neotoma_batches.py
  02_fix_missing_coordinates.py
  03_nordic_analysis\
    nordic_pollen_data.xlsx
  03_nordic_analysis.py
  04_exploratory_data_analysis\
    data_quality_report.xlsx
    EDA_SUMMARY_REPORT.txt
    site_statistics.xlsx
    spatial_basemap_ALL_IN_ONE_VERTICAL.tif
    spatial_basemap_diversity.tif
    spatial_basemap_samples.tif
    spatial_basemap_temporal.tif
    temporal_full_age_distribution.tif
  04_exploratory_data_analysis.py
  05_biotic_time_series\
    per_sample_metrics.csv
    plot_openness_period_loess.png
    plot_rate_of_change_period.png
    plot_shannon_period_rolling.png
    regional_bins_period_subbins.csv
    regional_with_loess.csv
    rolling_window_medians.csv
    site_bin_metrics_period_subbins.csv
    summary.json
  05_biotic_time_series.py
  merged_pollen_records.txt
  merged_pollen_records_fixed.txt
  run.sh
```

### Workflow overview

- 00 – Download: Download Neotoma datasets and extract geospatial metadata.
- 01 – Merge: Merge downloaded batches into unified records.
- 02 – Fix: Resolve missing or inconsistent coordinates.
- 03 – Nordic analysis: Subset/transform for Nordic region specifics and export consolidated Excel.
- 04 – EDA: Generate data quality reports, statistics, temporal distributions, and spatial rasters/maps.
- 05 – Biotic time series: Compute per-sample and regional time series metrics, summarize, and plot.

### Getting started

1) Requirements
- Python 3.9+ recommended
- Packages: pandas, numpy, requests, openpyxl, scipy, matplotlib, seaborn, geopandas, rasterio, shapely, scikit-learn, statsmodels (exact set may vary per script)

2) Installation (virtual environment recommended)
```bash
python -m venv .venv
. .venv/Scripts/Activate.ps1   # PowerShell on Windows
pip install --upgrade pip
pip install -r requirements.txt
```

If `requirements.txt` is not present, install packages on-demand based on script errors, or export your current environment using:
```bash
pip freeze > requirements.txt
```

### Running the pipeline

Option A – One-shot script
```bash
bash run.sh
```

Option B – Step-by-step
```bash
python 00_neotoma_download_extract_geo.py
python 01_merge_neotoma_batches.py
python 02_fix_missing_coordinates.py
python 03_nordic_analysis.py
python 04_exploratory_data_analysis.py
python 05_biotic_time_series.py
```

Outputs include both tabular summaries (CSV/XLSX) and figures/rasters (PNG/TIF). Some intermediate files are kept in the same directories for transparency and reproducibility.

### Data sources and citation

- Neotoma Paleoecology Database: please cite the database and any constituent datasets following their guidelines.

### Reproducibility notes

- Scripts assume Windows paths as shown above; adjust paths if running on other OSes.
- If large raster or figure generation is slow, ensure GDAL/RasterIO and GEOS are installed correctly in your environment.
- If you modify configuration (e.g., geographic bounds, time bins), commit those changes with a clear message so results can be traced.

### Contributing

Pull requests are welcome. For significant changes, please open an issue to discuss what you would like to change and ensure consistency across the workflow steps.

### License

Add your preferred license (e.g., MIT) in a `LICENSE` file.




