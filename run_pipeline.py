#!/usr/bin/env python3
"""
Master Pipeline Runner for Cross-Disease Transcriptomic Meta-Analysis of Spondyloarthritis.

Runs all analysis phases sequentially or individual phases via command-line arguments.

Usage:
    python run_pipeline.py --all                 # Run entire pipeline
    python run_pipeline.py --phase 3             # Run phase 3 only
    python run_pipeline.py --phase 1 1b 2        # Run phases 1, 1b, and 2
    python run_pipeline.py --list                # List all available phases
    python run_pipeline.py --from-phase 4        # Run from phase 4 onwards

Author: Pawel Kuklik
Spondyloarthritis AI Computational Biology Institute, Lisbon, Portugal
"""

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# Pipeline root
PIPELINE_DIR = Path(__file__).parent.resolve()
SCRIPTS_DIR = PIPELINE_DIR / "scripts"
LOGS_DIR = PIPELINE_DIR / "logs"

# Ordered pipeline phases with their scripts and descriptions
PHASES = [
    ("1",             "phase1_dataset_discovery.py",          "Dataset Discovery (NCBI GEO query)"),
    ("1b",            "phase1b_curate_datasets.py",           "Dataset Curation (21 → 7 datasets)"),
    ("2",             "phase2_download_preprocess.py",        "Download & Preprocessing"),
    ("2b",            "phase2b_fix_groups_download_rnaseq.py","Fix Groups & Download RNA-seq"),
    ("2c",            "phase2c_process_rnaseq.py",            "Process RNA-seq Data"),
    ("3",             "phase3_differential_expression.py",    "Differential Expression (Welch's t-test)"),
    ("3b",            "phase3b_probe_mapping.py",             "Probe-to-Gene Mapping"),
    ("3b_quick",      "phase3b_quick_mapping.py",             "Quick Probe Mapping"),
    ("4",             "phase4_meta_analysis.py",              "Fisher's Meta-Analysis"),
    ("5",             "phase5_enrichment_and_figures.py",     "Functional Enrichment & Figures"),
    ("6",             "phase6_ppi_network.py",                "PPI Network Analysis (STRING)"),
    ("7",             "phase7_wgcna.py",                      "WGCNA Co-expression Network"),
    ("8",             "phase8_immune_deconvolution.py",       "Immune Cell Deconvolution (ssGSEA)"),
    ("8b",            "phase8b_fix_psa_deconvolution.py",     "Fix PsA Deconvolution"),
    ("8c",            "phase8c_update_crossdisease.py",       "Update Cross-Disease Comparisons"),
    ("sensitivity",   "sensitivity_analysis.py",              "Sensitivity Analysis"),
    ("supplementary", "generate_supplementary_tables.py",     "Generate Supplementary Tables"),
    ("abstract",      "graphical_abstract.py",                "Graphical Abstract"),
    ("figures",       "inject_figures.py",                    "Inject Figures into Manuscript"),
    ("figures_v2",    "inject_figures_v2.py",                 "Inject Figures v2 (Revised)"),
]


def get_phase_ids():
    """Return list of valid phase identifiers."""
    return [p[0] for p in PHASES]


def get_phase_by_id(phase_id):
    """Look up a phase by its identifier."""
    for pid, script, desc in PHASES:
        if pid == phase_id:
            return pid, script, desc
    return None


def format_duration(seconds):
    """Format elapsed time as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    secs = seconds % 60
    if minutes < 60:
        return f"{minutes}m {secs:.0f}s"
    hours = int(minutes // 60)
    mins = minutes % 60
    return f"{hours}h {mins}m {secs:.0f}s"


def print_header(phase_id, description, phase_num, total):
    """Print a formatted phase header."""
    width = 70
    print()
    print("=" * width)
    print(f"  Phase {phase_id}: {description}")
    print(f"  [{phase_num}/{total}] | Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * width)
    print()


def print_summary(results):
    """Print pipeline execution summary."""
    width = 70
    print()
    print("=" * width)
    print("  PIPELINE EXECUTION SUMMARY")
    print("=" * width)
    print()
    print(f"  {'Phase':<15} {'Status':<12} {'Duration':<12} Description")
    print(f"  {'-'*13:<15} {'-'*10:<12} {'-'*10:<12} {'-'*30}")

    total_time = 0
    passed = 0
    failed = 0

    for phase_id, desc, status, duration in results:
        status_str = "PASSED" if status == 0 else f"FAILED ({status})"
        dur_str = format_duration(duration)
        total_time += duration

        if status == 0:
            passed += 1
        else:
            failed += 1

        print(f"  {phase_id:<15} {status_str:<12} {dur_str:<12} {desc}")

    print()
    print(f"  Total: {passed} passed, {failed} failed, "
          f"{format_duration(total_time)} elapsed")
    print("=" * width)
    print()


def run_phase(phase_id, script, description, phase_num, total, log_dir=None):
    """Run a single pipeline phase script and return (status_code, duration)."""
    print_header(phase_id, description, phase_num, total)

    script_path = SCRIPTS_DIR / script
    if not script_path.exists():
        print(f"  ERROR: Script not found: {script_path}")
        return 1, 0.0

    cmd = [sys.executable, str(script_path)]
    start = time.time()

    log_file = None
    if log_dir:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"phase_{phase_id}.log"
        log_file = open(log_path, "w")

    try:
        result = subprocess.run(
            cmd,
            cwd=str(PIPELINE_DIR),
            stdout=log_file if log_file else None,
            stderr=subprocess.STDOUT if log_file else None,
            timeout=7200,  # 2 hour timeout per phase
        )
        status = result.returncode
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT: Phase {phase_id} exceeded 2-hour limit")
        status = -1
    except Exception as e:
        print(f"  ERROR running phase {phase_id}: {e}")
        status = -2
    finally:
        if log_file:
            log_file.close()

    duration = time.time() - start
    status_msg = "completed successfully" if status == 0 else f"failed (exit code {status})"
    print(f"\n  Phase {phase_id} {status_msg} in {format_duration(duration)}")

    if log_file:
        print(f"  Log: {log_path}")

    return status, duration


def list_phases():
    """Print all available pipeline phases."""
    print("\nAvailable Pipeline Phases:")
    print(f"  {'ID':<15} {'Script':<42} Description")
    print(f"  {'-'*13:<15} {'-'*40:<42} {'-'*30}")
    for phase_id, script, desc in PHASES:
        print(f"  {phase_id:<15} {script:<42} {desc}")
    print(f"\nTotal: {len(PHASES)} phases")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Cross-Disease Transcriptomic Meta-Analysis Pipeline Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py --all                  Run entire pipeline
  python run_pipeline.py --phase 3              Run phase 3 only
  python run_pipeline.py --phase 1 1b 2 3       Run phases 1, 1b, 2, 3
  python run_pipeline.py --from-phase 4         Run from phase 4 to end
  python run_pipeline.py --list                 List all phases
  python run_pipeline.py --all --log            Save output to log files
        """,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--all",
        action="store_true",
        help="Run all pipeline phases sequentially",
    )
    group.add_argument(
        "--phase",
        nargs="+",
        metavar="ID",
        help="Run specific phase(s) by ID (e.g., 1, 3b, sensitivity)",
    )
    group.add_argument(
        "--from-phase",
        metavar="ID",
        help="Run from specified phase through end of pipeline",
    )
    group.add_argument(
        "--list",
        action="store_true",
        help="List all available pipeline phases",
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Save phase output to log files in logs/ directory",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop pipeline on first phase failure (default: continue)",
    )

    args = parser.parse_args()

    if args.list:
        list_phases()
        return

    # Determine which phases to run
    if args.all:
        phases_to_run = PHASES[:]
    elif args.phase:
        phases_to_run = []
        for pid in args.phase:
            phase = get_phase_by_id(pid)
            if phase is None:
                print(f"Error: Unknown phase '{pid}'. Use --list to see available phases.")
                sys.exit(1)
            phases_to_run.append(phase)
    elif args.from_phase:
        phase_ids = get_phase_ids()
        if args.from_phase not in phase_ids:
            print(f"Error: Unknown phase '{args.from_phase}'. Use --list to see available phases.")
            sys.exit(1)
        start_idx = phase_ids.index(args.from_phase)
        phases_to_run = PHASES[start_idx:]

    # Log directory
    log_dir = LOGS_DIR if args.log else None

    # Pipeline banner
    print()
    print("*" * 70)
    print("  Cross-Disease Transcriptomic Meta-Analysis of Spondyloarthritis")
    print("  Pipeline Runner")
    print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Phases to run: {len(phases_to_run)}")
    print("*" * 70)

    # Run phases
    results = []
    total = len(phases_to_run)

    for i, (phase_id, script, desc) in enumerate(phases_to_run, 1):
        status, duration = run_phase(phase_id, script, desc, i, total, log_dir)
        results.append((phase_id, desc, status, duration))

        if status != 0 and args.stop_on_error:
            print(f"\n  Stopping pipeline due to error in phase {phase_id}")
            break

    # Summary
    print_summary(results)

    # Exit with error if any phase failed
    failed = sum(1 for _, _, s, _ in results if s != 0)
    if failed:
        print(f"WARNING: {failed} phase(s) failed. Check logs for details.")
        sys.exit(1)

    print("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
