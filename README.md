## Tundra To Forest

Workflow to download, clean, analyze, and visualize pollen records for tundra→forest transitions.

### Quick start

1) Create a virtual environment and install deps
```bash
python -m venv .venv
. .venv/Scripts/Activate.ps1
pip install --upgrade pip
pip install -r requirements.txt
```

2) Run the pipeline (step-by-step)
```bash
python 00_neotoma_download_extract_geo.py
python 01_merge_neotoma_batches.py
python 02_fix_missing_coordinates.py
python 03_nordic_analysis.py
python 04_exploratory_data_analysis.py
python 05_biotic_time_series.py
```

Optional one-shot
```bash
bash run.sh
```

### Outputs
- Data summaries (CSV/XLSX)
- Figures and rasters (PNG/TIF)

### Notes
- Large raster/image assets are tracked with Git LFS.
- Some large intermediate text files may be excluded; regenerate them by re-running the steps above.

### Citation
Please cite Neotoma Paleoecology Database and any constituent datasets used.

### License
Add a `LICENSE` file (e.g., MIT) if you plan to distribute.