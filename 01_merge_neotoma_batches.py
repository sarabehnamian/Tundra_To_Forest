#!/usr/bin/env python3
"""
merge_neotoma_batches.py
Merge all pollen_records.txt files from batch subdirectories into a single consolidated file.
"""

import pandas as pd
from pathlib import Path
import argparse
import sys

def merge_batch_results(batch_dir: Path, output_file: Path, debug: bool = False) -> int:
    """
    Merge all pollen_records.txt files from batch subdirectories.
    
    Args:
        batch_dir: Root directory containing batch_* subdirectories
        output_file: Path to write the merged output
        debug: Enable verbose output
        
    Returns:
        Total number of rows in merged dataset
    """
    
    # Find all batch subdirectories
    batch_dirs = sorted([d for d in batch_dir.iterdir() if d.is_dir() and d.name.startswith("batch_")])
    
    if not batch_dirs:
        print(f"ERROR: No batch_* directories found in {batch_dir}")
        return 0
    
    print(f"Found {len(batch_dirs)} batch directories to process")
    
    # Collect all dataframes
    all_dfs = []
    total_rows = 0
    missing_files = []
    empty_files = []
    
    for i, bdir in enumerate(batch_dirs, 1):
        pollen_file = bdir / "pollen_records.txt"
        
        # Check if file exists
        if not pollen_file.exists():
            missing_files.append(bdir.name)
            if debug:
                print(f"[{i}/{len(batch_dirs)}] {bdir.name}: MISSING pollen_records.txt")
            continue
        
        # Check if DONE marker exists (optional check)
        done_markers = list(bdir.glob("DONE_*"))
        has_done = len(done_markers) > 0
        
        try:
            # Read the CSV file
            df = pd.read_csv(pollen_file, low_memory=False)
            
            if len(df) == 0:
                empty_files.append(bdir.name)
                if debug:
                    print(f"[{i}/{len(batch_dirs)}] {bdir.name}: EMPTY (0 rows)")
                continue
            
            all_dfs.append(df)
            total_rows += len(df)
            
            status = "COMPLETE" if has_done else "NO_DONE_MARKER"
            print(f"[{i}/{len(batch_dirs)}] {bdir.name}: {len(df):,} rows ({status})")
            
        except Exception as e:
            print(f"[{i}/{len(batch_dirs)}] {bdir.name}: ERROR reading file - {e}")
            continue
    
    # Report any issues
    if missing_files:
        print(f"\nWARNING: {len(missing_files)} batches missing pollen_records.txt:")
        for mf in missing_files[:10]:  # Show first 10
            print(f"  - {mf}")
        if len(missing_files) > 10:
            print(f"  ... and {len(missing_files) - 10} more")
    
    if empty_files:
        print(f"\nINFO: {len(empty_files)} batches had empty pollen_records.txt files")
    
    # Merge all dataframes
    if not all_dfs:
        print("\nERROR: No valid data files found to merge!")
        return 0
    
    print(f"\nMerging {len(all_dfs)} dataframes with {total_rows:,} total rows...")
    merged_df = pd.concat(all_dfs, ignore_index=True)
    
    # Check for duplicates based on key columns (optional)
    if debug:
        key_cols = ['dataset_id', 'sample_id', 'taxon']
        available_keys = [col for col in key_cols if col in merged_df.columns]
        if available_keys:
            n_duplicates = merged_df.duplicated(subset=available_keys).sum()
            if n_duplicates > 0:
                print(f"WARNING: Found {n_duplicates:,} potential duplicate rows based on {available_keys}")
    
    # Sort by dataset_id and sample_id for consistency (optional)
    if 'dataset_id' in merged_df.columns:
        sort_cols = ['dataset_id']
        if 'sample_id' in merged_df.columns:
            sort_cols.append('sample_id')
        merged_df = merged_df.sort_values(sort_cols)
        print(f"Sorted by: {sort_cols}")
    
    # Write the merged file
    print(f"\nWriting merged data to: {output_file}")
    merged_df.to_csv(output_file, index=False)
    
    # Final statistics
    print(f"\n{'='*60}")
    print(f"MERGE COMPLETE")
    print(f"{'='*60}")
    print(f"Total batches processed: {len(all_dfs)}/{len(batch_dirs)}")
    print(f"Total rows merged: {len(merged_df):,}")
    print(f"Output file: {output_file}")
    print(f"Output file size: {output_file.stat().st_size / (1024**2):.1f} MB")
    
    # Column summary
    if debug:
        print(f"\nColumns in merged dataset ({len(merged_df.columns)}):")
        for col in merged_df.columns:
            non_null = merged_df[col].notna().sum()
            pct = 100 * non_null / len(merged_df)
            print(f"  - {col}: {non_null:,} non-null ({pct:.1f}%)")
    
    return len(merged_df)


def verify_merge(merged_file: Path, batch_dir: Path) -> None:
    """
    Verify that the merge was successful by comparing row counts.
    """
    print("\n" + "="*60)
    print("VERIFICATION")
    print("="*60)
    
    # Count total rows in individual files
    expected_total = 0
    for bdir in batch_dir.iterdir():
        if not bdir.is_dir() or not bdir.name.startswith("batch_"):
            continue
        pollen_file = bdir / "pollen_records.txt"
        if pollen_file.exists():
            try:
                df = pd.read_csv(pollen_file, nrows=0)  # Just get column count
                row_count = sum(1 for _ in open(pollen_file)) - 1  # Subtract header
                expected_total += row_count
            except:
                pass
    
    # Count rows in merged file
    merged_count = sum(1 for _ in open(merged_file)) - 1  # Subtract header
    
    print(f"Expected total rows (sum of batches): {expected_total:,}")
    print(f"Actual rows in merged file: {merged_count:,}")
    
    if expected_total == merged_count:
        print("✓ Verification PASSED: Row counts match!")
    else:
        diff = merged_count - expected_total
        print(f"⚠ Verification WARNING: Difference of {diff:+,} rows")
        if diff < 0:
            print("  (Some rows may have been lost during merge)")
        else:
            print("  (This could be due to header rows or formatting)")


def main():
    parser = argparse.ArgumentParser(
        description="Merge all pollen_records.txt files from Neotoma batch directories"
    )
    parser.add_argument(
        "--batch-dir",
        type=Path,
        default=Path("neotoma_out_batches"),
        help="Root directory containing batch_* subdirectories (default: neotoma_out_batches)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("merged_pollen_records.txt"),
        help="Output file path (default: merged_pollen_records.txt)"
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify the merge by comparing row counts"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable verbose debug output"
    )
    
    args = parser.parse_args()
    
    # Check if batch directory exists
    if not args.batch_dir.exists():
        print(f"ERROR: Batch directory not found: {args.batch_dir}")
        sys.exit(1)
    
    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # Perform the merge
    n_rows = merge_batch_results(args.batch_dir, args.output, debug=args.debug)
    
    if n_rows > 0 and args.verify:
        verify_merge(args.output, args.batch_dir)
    
    sys.exit(0 if n_rows > 0 else 1)


if __name__ == "__main__":
    main()