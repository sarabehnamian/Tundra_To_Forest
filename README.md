## From tundra to forest and back: 21,000 years of Nordic landscape transformation

Pipeline to download, clean, analyze, and visualize pollen records for reconstructing Nordic vegetation dynamics over the last 21,000 years.

### Quick start

1) Create a virtual environment and install dependencies
```bash
python -m venv .venv
. .venv/Scripts/Activate.ps1
pip install --upgrade pip
pip install -r requirements.txt
```

2) Run the pipeline (step-by-step)
```bash
python 00_neotoma_download_extract_geo.py         # Download from Neotoma and extract geospatial metadata
python 01_merge_neotoma_batches.py                # Merge downloaded batches into unified records
python 02_fix_missing_coordinates.py              # Resolve missing/inconsistent coordinates
python 03_nordic_analysis.py                      # Nordic-specific subsetting and transforms
python 04_exploratory_data_analysis.py            # Data QC, EDA summaries, site statistics, temporal distributions
python 05_biotic_time_series.py                   # Openness, diversity (H', richness, evenness), turnover time series
```

Optional one-shot
```bash
bash run.sh
```

### What each script does
- `00_neotoma_download_extract_geo.py`: Downloads pollen datasets from Neotoma, extracts site-level geospatial attributes, and writes standardized tables for downstream steps.
- `01_merge_neotoma_batches.py`: Consolidates multiple Neotoma pulls into a single record table, harmonizing field names and types.
- `02_fix_missing_coordinates.py`: Detects and fixes missing or inconsistent latitude/longitude entries where possible.
- `03_nordic_analysis.py`: Applies Nordic domain filters and transformations; prepares region-focused aggregates and exports consolidated Excel where needed.
- `04_exploratory_data_analysis.py`: Produces quality-control summaries, site statistics, temporal age distributions, and other EDA artifacts.
- `05_biotic_time_series.py`: Computes ecological indicators (openness = NAP/(AP+NAP), diversity metrics, chord-distance turnover), applies adaptive 500-year bins, bootstrap uncertainty, and optional LOESS smoothing.

### Outputs
- Tabular outputs (CSV/XLSX) for summaries and intermediate aggregates
- Figures (PNG) and, when enabled, rasters (TIF) for spatial/temporal visualizations

Note: Very large intermediate artifacts (e.g., giant TXT, TIF rasters) may be excluded from version control. Recreate them by re-running the corresponding steps.

### How to cite
From tundra to forest and back: 21,000 years of Nordic landscape transformation

Authors: Sara Behnamian¹, Fatemeh Fogh²

Affiliations: ¹Globe Institute, University of Copenhagen, Øster Voldgade 5–7, 1350 Copenhagen K, Denmark ²Department of Mathematical Sciences, Florida Atlantic University, Boca Raton, FL 33431, USA

Corresponding Author: Sara Behnamian  Email: sara.behnamian@sund.ku.dk

Abstract

We present a quantitative reconstruction of 21,000 years of Nordic vegetation dynamics using 271,462 pollen records from 294 sites compiled from the Neotoma Paleoecology Database. After systematic quality control, taxonomic filtering, and harmonized chronologies, we quantified three key ecological metrics: landscape openness (non-arboreal pollen ratio), taxonomic diversity (richness, Shannon’s H′, evenness), and vegetation turnover (chord distance). Adaptive 500-year temporal binning, bootstrap uncertainty (1,000 iterations), and LOESS smoothing (span = 0.25) were used to resolve long-term trends.

Results show a transition from open Late Pleistocene tundra–steppe landscapes (median openness = 0.268) to densely forested Mid-Holocene conditions (0.039), representing 85 % landscape closure during post-glacial forest expansion. Diversity and compositional turnover peaked during deglaciation, declined under mid-Holocene forest homogenization, and rose again through the Late Holocene. Historical re-opening (openness = 0.114) indicates widespread anthropogenic deforestation.

This synthesis provides a reproducible, data-driven framework linking large-scale vegetation trajectories to climate and land-use change. By harmonizing regional pollen archives through standardized processing, it establishes a robust baseline for evaluating future Nordic ecosystem responses under accelerating environmental pressure.

### Data sources
- Neotoma Paleoecology Database (cite according to their guidelines)

### License
MIT License - see LICENSE file for details.