#!/usr/bin/env python3
"""
02_fix_missing_coordinates.py
Fix missing coordinates AND generate detailed reports, saving everything in a subfolder.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from datetime import datetime

def find_pollen_file():
    """
    Automatically find the pollen records file to process.
    """
    search_patterns = [
        "merged_pollen_records.txt",
        "pollen_records.txt",
        "*pollen*.txt",
        "*pollen*.csv"
    ]
    
    current_dir = Path(".")
    
    for pattern in search_patterns:
        files = list(current_dir.glob(pattern))
        # Exclude files that end with _fixed or contain 'missing'
        files = [f for f in files if not f.stem.endswith('_fixed') and 'missing' not in f.stem.lower()]
        if files:
            if len(files) > 1:
                files.sort(key=lambda x: x.stat().st_size, reverse=True)
            return files[0]
    
    print("ERROR: Could not find pollen records file!")
    sys.exit(1)

def main():
    """
    Main function - fix coordinates and generate reports.
    """
    
    # Find the input file
    input_file = find_pollen_file()
    print(f"Found file: {input_file}")
    print(f"File size: {input_file.stat().st_size / (1024**2):.1f} MB")
    
    # Create output folder
    output_folder = Path("02_fix_missing_coordinates")
    output_folder.mkdir(exist_ok=True)
    print(f"\nCreated output folder: {output_folder}")
    
    # Read the data
    print(f"\nReading data...")
    df = pd.read_csv(input_file, low_memory=False)
    print(f"Total records loaded: {len(df):,}")
    
    # Identify missing coordinates
    lat_missing = df['latitude'].isna() | (df['latitude'] == '') | (df['latitude'].astype(str).str.strip() == '')
    lon_missing = df['longitude'].isna() | (df['longitude'] == '') | (df['longitude'].astype(str).str.strip() == '')
    coords_missing = lat_missing | lon_missing
    
    n_missing = coords_missing.sum()
    print(f"\nRecords with missing coordinates: {n_missing:,} ({100*n_missing/len(df):.2f}%)")
    
    # Save all records with missing coordinates BEFORE fixing
    df_missing_original = df[coords_missing].copy()
    missing_records_file = output_folder / "records_with_missing_coordinates.txt"
    print(f"\nSaving {len(df_missing_original):,} records with missing coordinates...")
    df_missing_original.to_csv(missing_records_file, index=False, sep='\t')
    print(f"Saved to: {missing_records_file}")
    
    if n_missing == 0:
        print("No missing coordinates found. Nothing to fix!")
        return
    
    # Get unique sites with missing coordinates
    sites_with_missing = df[coords_missing]['site_name'].unique()
    sites_with_missing = [s for s in sites_with_missing if pd.notna(s)]
    print(f"\nUnique sites with missing coordinates: {len(sites_with_missing)}")
    
    # Build lookup table from records with coordinates
    print("\nBuilding coordinate lookup table...")
    has_coords = ~coords_missing & df['latitude'].notna() & df['longitude'].notna()
    
    coords_lookup = {}
    sites_with_conflicts = []
    sites_no_coords = []
    
    for site in sites_with_missing:
        site_mask = (df['site_name'] == site) & has_coords
        site_data = df[site_mask]
        
        if len(site_data) == 0:
            sites_no_coords.append(site)
            continue
        
        site_data_rounded = site_data.copy()
        site_data_rounded['lat_round'] = site_data_rounded['latitude'].round(5)
        site_data_rounded['lon_round'] = site_data_rounded['longitude'].round(5)
        
        unique_coords = site_data_rounded[['lat_round', 'lon_round']].drop_duplicates()
        
        if len(unique_coords) == 1:
            coords_lookup[site] = (unique_coords.iloc[0]['lat_round'], unique_coords.iloc[0]['lon_round'])
        elif len(unique_coords) > 1:
            sites_with_conflicts.append(site)
    
    print(f"Found reliable coordinates for {len(coords_lookup)} sites")
    print(f"Sites with conflicting coordinates (skipped): {len(sites_with_conflicts)}")
    print(f"Sites with no coordinate data available: {len(sites_no_coords)}")
    
    # Apply the fixes
    print("\nApplying coordinate fixes...")
    fixed_count = 0
    df_fixed = df.copy()
    fixed_records = []
    
    for idx in df[coords_missing].index:
        site = df.loc[idx, 'site_name']
        
        if site in coords_lookup:
            lat, lon = coords_lookup[site]
            df_fixed.loc[idx, 'latitude'] = lat
            df_fixed.loc[idx, 'longitude'] = lon
            fixed_records.append(idx)
            fixed_count += 1
        
        if fixed_count > 0 and fixed_count % 100000 == 0:
            print(f"  Progress: {fixed_count:,} records fixed...")
    
    print(f"Fixed coordinates for {fixed_count:,} records")
    
    # Check remaining missing
    lat_missing_after = df_fixed['latitude'].isna()
    lon_missing_after = df_fixed['longitude'].isna()
    still_missing = lat_missing_after | lon_missing_after
    n_still_missing = still_missing.sum()
    
    print(f"Records still missing coordinates: {n_still_missing:,}")
    
    # Save the fixed complete dataset
    output_file = output_folder / f"{input_file.stem}_fixed.txt"
    print(f"\nSaving fixed dataset to: {output_file}")
    df_fixed.to_csv(output_file, index=False, sep='\t')
    
    # Save records that were fixed
    if fixed_records:
        df_fixed_records = df_fixed.loc[fixed_records]
        fixed_records_file = output_folder / "records_that_were_fixed.txt"
        df_fixed_records.to_csv(fixed_records_file, index=False, sep='\t')
        print(f"Saved fixed records to: {fixed_records_file}")
    
    # Save records still missing coordinates
    df_still_missing = df_fixed[still_missing]
    still_missing_file = output_folder / "records_still_missing_coordinates.txt"
    df_still_missing.to_csv(still_missing_file, index=False, sep='\t')
    print(f"Saved still missing records to: {still_missing_file}")
    
    # Generate detailed reports
    print("\nGenerating reports...")
    
    # Site-level analysis
    sites_summary = df_missing_original.groupby('site_name').agg({
        'dataset_id': ['count', 'nunique', lambda x: list(x.unique())[:10]],
        'sample_id': 'nunique',
        'taxon': 'nunique'
    }).round(2)
    sites_summary.columns = ['record_count', 'n_datasets', 'dataset_ids_sample', 'n_samples', 'n_taxa']
    sites_summary = sites_summary.sort_values('record_count', ascending=False)
    
    # Dataset-level analysis
    dataset_summary = df_missing_original.groupby('dataset_id').agg({
        'site_name': 'first',
        'latitude': 'count',
        'sample_id': 'nunique',
        'taxon': 'nunique'
    })
    dataset_summary.columns = ['site_name', 'missing_records', 'n_samples', 'n_taxa']
    dataset_summary = dataset_summary.sort_values('missing_records', ascending=False)
    
    # Create comprehensive text report
    report_file = output_folder / "REPORT_missing_coordinates.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("MISSING COORDINATES ANALYSIS AND FIX REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n\n")
        
        f.write("SUMMARY\n")
        f.write("-"*40 + "\n")
        f.write(f"Input file: {input_file.name}\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Records with missing coordinates: {n_missing:,} ({100*n_missing/len(df):.2f}%)\n")
        f.write(f"Records fixed: {fixed_count:,}\n")
        f.write(f"Records still missing: {n_still_missing:,}\n")
        f.write(f"Fix success rate: {100*fixed_count/n_missing:.1f}%\n\n")
        
        f.write("SITE ANALYSIS\n")
        f.write("-"*40 + "\n")
        f.write(f"Total sites with missing coordinates: {len(sites_with_missing)}\n")
        f.write(f"Sites with coordinates found: {len(coords_lookup)}\n")
        f.write(f"Sites with conflicting coordinates: {len(sites_with_conflicts)}\n")
        f.write(f"Sites with no coordinate data: {len(sites_no_coords)}\n\n")
        
        f.write("TOP 30 SITES WITH MISSING COORDINATES\n")
        f.write("-"*40 + "\n")
        for idx, (site_name, row) in enumerate(sites_summary.head(30).iterrows(), 1):
            status = "FIXED" if site_name in coords_lookup else "CONFLICT" if site_name in sites_with_conflicts else "NO DATA"
            f.write(f"{idx:3}. {site_name:<50} {row['record_count']:8.0f} records  [{status}]\n")
        
        f.write("\n\nTOP 30 DATASETS WITH MISSING COORDINATES\n")
        f.write("-"*40 + "\n")
        for idx, (dataset_id, row) in enumerate(dataset_summary.head(30).iterrows(), 1):
            f.write(f"{idx:3}. Dataset {dataset_id:6}: {row['missing_records']:8} records  (Site: {row['site_name']})\n")
        
        if sites_with_conflicts:
            f.write("\n\nSITES WITH CONFLICTING COORDINATES (NOT FIXED)\n")
            f.write("-"*40 + "\n")
            for idx, site in enumerate(sites_with_conflicts[:50], 1):
                f.write(f"{idx:3}. {site}\n")
    
    print(f"Saved report to: {report_file}")
    
    # Save site summary as tab-delimited
    site_report_file = output_folder / "sites_missing_coordinates_details.txt"
    sites_summary.to_csv(site_report_file, sep='\t')
    print(f"Saved site details to: {site_report_file}")
    
    # Save dataset summary as tab-delimited
    dataset_report_file = output_folder / "datasets_missing_coordinates_details.txt"
    dataset_summary.to_csv(dataset_report_file, sep='\t')
    print(f"Saved dataset details to: {dataset_report_file}")
    
    # Print final summary
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print(f"All output saved to folder: {output_folder}/")
    print("\nGenerated files:")
    print(f"1. {output_file.name} - Complete dataset with fixes applied")
    print(f"2. records_with_missing_coordinates.txt - Original records missing coords")
    print(f"3. records_that_were_fixed.txt - Records that got fixed")
    print(f"4. records_still_missing_coordinates.txt - Records still missing coords")
    print(f"5. REPORT_missing_coordinates.txt - Main analysis report")
    print(f"6. sites_missing_coordinates_details.txt - Detailed site analysis")
    print(f"7. datasets_missing_coordinates_details.txt - Detailed dataset analysis")
    print("\n" + "="*80)

if __name__ == "__main__":
    main()