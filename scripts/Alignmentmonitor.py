#!/usr/bin/env python3

import os
import time
import csv
from pathlib import Path
from datetime import datetime

# Config
TOTAL_SAMPLES = 734
OUTPUT_BASE = "/path/to/your/results"

# Define how to detect a completed sample for each aligner.
# For each aligner, provide:
#   - 'dir': subdirectory under OUTPUT_BASE where results live
#   - 'marker': a file or directory name that indicates completion
#               (e.g., 'Aligned.sortedByCoord.out.bam', 'quant.sf', or a '.done' file)
ALIGNERS = {
    'STAR_with_Load': {
        'dir': 'STAR_genomeLoad',
        'marker': 'Aligned.sortedByCoord.out.bam'  
    },
    'Salmon': {
        'dir': 'Salmon_quants',
        'marker': 'quant.sf'                         
    },
    'STAR_no_Load': {
        'dir': 'STAR_noLoad',
        'marker': 'Aligned.sortedByCoord.out.bam'
    },
    'TopHat2': {
        'dir': 'TopHat2_results',
        'marker': 'accepted_hits.bam'                
    },
}

# Interval between checks = 1800 seconds (30 minutes)
CHECK_INTERVAL = 1800

# Output CSV file
OUTPUT_CSV = 'aligner_progress_real.csv'

# Recomennded to test before editing (or not to edit at all), can cause unintended breaks in driver (Known issue, cause is unknown)

def count_completed_samples(base_dir, sub_dir, marker):
    """Count how many sample directories contain the marker file."""
    target_dir = Path(base_dir) / sub_dir
    if not target_dir.exists():
        print(f"Warning: Directory {target_dir} does not exist yet.")
        return 0

    count = 0
    # Assumes each sample has its own subdirectory
    for sample_dir in target_dir.iterdir():
        if sample_dir.is_dir():
            marker_path = sample_dir / marker
            if marker_path.exists():
                count += 1
    return count

def main():
    print(f"Starting aligner progress monitor.")
    print(f"Checking every {CHECK_INTERVAL} seconds ({CHECK_INTERVAL/60:.1f} minutes).")
    print(f"Press Ctrl+C to stop.\n")

    # Initialize CSV file
    header = ['Timestamp', 'Elapsed_Hours'] + list(ALIGNERS.keys())
    with open(OUTPUT_CSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

    start_time = time.time()

    try:
        while True:
            elapsed_sec = time.time() - start_time
            elapsed_hours = elapsed_sec / 3600.0
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

            # Count completions for each aligner
            counts = {}
            all_done = True
            for name, cfg in ALIGNERS.items():
                cnt = count_completed_samples(OUTPUT_BASE, cfg['dir'], cfg['marker'])
                counts[name] = cnt
                if cnt < TOTAL_SAMPLES:
                    all_done = False

            # Write to CSV
            row = [timestamp, f"{elapsed_hours:.2f}"] + [counts[name] for name in ALIGNERS]
            with open(OUTPUT_CSV, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(row)

            # Print status to console
            status = ' | '.join([f"{name}: {counts[name]}/{TOTAL_SAMPLES}" for name in ALIGNERS])
            print(f"[{timestamp}] Elapsed: {elapsed_hours:.2f}h | {status}")

            # Exit if all samples are done for all aligners
            if all_done:
                print("\nAll aligners have finished processing all samples!")
                break

            # Wait for next interval
            time.sleep(CHECK_INTERVAL)

    except KeyboardInterrupt:
        print("\nMonitoring stopped by user.")
        print(f"Partial progress saved to {OUTPUT_CSV}")

if __name__ == '__main__':
    main()