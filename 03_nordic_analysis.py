#!/usr/bin/env python3
"""
03_nordic_analysis.py
Filter pollen data for Nordic countries/territories (with valid age data) and save as Excel.

Included regions:
- Denmark (incl. Bornholm)
- Finland (incl. Åland)
- Iceland
- Norway (mainland, Svalbard, Jan Mayen)
- Sweden
- Faroe Islands (DK)

Input : merged_pollen_records_fixed.txt (tab-separated)
Output: 03_nordic_analysis/nordic_pollen_data.xlsx
"""

import pandas as pd
from pathlib import Path
import sys


def is_in_nordic(lat, lon):
    """Return True if (lat, lon) lies in any Nordic country/territory box (excludes Greenland)."""
    if pd.isna(lat) or pd.isna(lon):
        return False

    # --- Core countries (bounding boxes) ---
    in_denmark = (54.4 <= lat <= 57.95) and (7.8 <= lon <= 15.2)   # Jutland/Funen/Zealand/Bornholm
    in_sweden  = (55.3 <= lat <= 69.1)  and (11.1 <= lon <= 24.2)
    in_norway  = (57.9 <= lat <= 71.3)  and (4.4  <= lon <= 31.2)  # mainland
    in_finland = (59.5 <= lat <= 70.2)  and (20.5 <= lon <= 31.6)
    in_iceland = (63.1 <= lat <= 66.7)  and (-24.6 <= lon <= -13.3)

    # --- Territories ---
    in_aland     = (59.9 <= lat <= 60.6) and (19.2 <= lon <= 21.3)
    in_faroe     = (61.3 <= lat <= 62.5) and (-7.9 <= lon <= -6.0)
    in_svalbard  = (74.0 <= lat <= 81.0) and (10.0 <= lon <= 35.0)
    in_jan_mayen = (70.5 <= lat <= 71.2) and (-9.5 <= lon <= -7.7)

    return any([
        in_denmark, in_sweden, in_norway, in_finland, in_iceland,
        in_aland, in_faroe, in_svalbard, in_jan_mayen
    ])


def main():
    print("\nFILTERING NORDIC POLLEN DATA")
    print("=" * 50)

    # Input file
    input_file = Path("merged_pollen_records_fixed.txt")
    if not input_file.exists():
        print(f"ERROR: Cannot find {input_file}")
        sys.exit(1)

    # Create output directory
    output_dir = Path("03_nordic_analysis")
    output_dir.mkdir(exist_ok=True)

    # Load data
    print(f"Loading data from: {input_file}")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    print(f"Total records loaded: {len(df):,}")

    # Required columns
    required_cols = ["latitude", "longitude", "age_BP", "element"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"ERROR: Missing required columns: {', '.join(missing)}")
        sys.exit(1)

    # Filter for Nordic region
    print("Filtering for Nordic coordinates (Greenland EXCLUDED)...")
    df["in_nordic"] = df.apply(
        lambda r: is_in_nordic(r["latitude"], r["longitude"]),
        axis=1
    )
    df_nordic = df[df["in_nordic"]].copy().drop(columns=["in_nordic"])
    print(f"Nordic records (before age filter): {len(df_nordic):,}")

    # Valid age_BP
    print("Filtering for records with valid age_BP...")
    df_nordic = df_nordic[df_nordic["age_BP"].notna()].copy()
    print(f"Nordic records with valid age: {len(df_nordic):,}")

    # Pollen-only (robust to NaN/mixed case)
    print("Filtering for pollen records only...")
    df_nordic["element_lower"] = df_nordic["element"].astype(str).str.strip().str.lower()
    df_nordic = df_nordic[df_nordic["element_lower"] == "pollen"].copy()
    df_nordic.drop(columns=["element_lower"], inplace=True)
    print(f"Nordic pollen records: {len(df_nordic):,}")
    print(f"Percentage of total: {100 * len(df_nordic) / max(len(df), 1):.2f}%")

    # Drop unwanted columns (only if present)
    columns_to_drop = [
        "date_BP", "age_older", "age_younger", "units",
        "year_CE", "year_older_CE", "year_younger_CE", "year_midpoint_CE",
        "chronology_name", "chronology_id", "element", "age_source",
        "collection_date"
    ]
    drop_existing = [c for c in columns_to_drop if c in df_nordic.columns]
    if drop_existing:
        print(f"Dropping columns: {', '.join(drop_existing)}")
        df_nordic.drop(columns=drop_existing, inplace=True)

    # Save as Excel (Excel row limit aware)
    excel_path = output_dir / "nordic_pollen_data.xlsx"
    print(f"\nSaving to Excel: {excel_path}")
    max_rows = 1_048_576
    if len(df_nordic) >= max_rows:
        print(f"WARNING: {len(df_nordic):,} rows exceed Excel limit ({max_rows:,}). Saving first {max_rows-1:,} rows.")
        df_nordic.head(max_rows - 1).to_excel(excel_path, index=False)
        saved_n = max_rows - 1
    else:
        df_nordic.to_excel(excel_path, index=False)
        saved_n = len(df_nordic)

    print(f"Done. Saved {saved_n:,} records to {excel_path}")


if __name__ == "__main__":
    main()
