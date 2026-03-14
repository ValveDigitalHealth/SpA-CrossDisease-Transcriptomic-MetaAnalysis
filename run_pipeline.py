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

PIPELINE_DIR = Path(__file__).parent.resolve()
SCRIPTS_DIR = PIPELINE_DIR / "scripts"
LOGS_DIR = PIPELINE_DIR / "logs"

PHASES = [
    ("1",             "phase1_dataset_discovery.py",          "Dataset Discovery (NCBI GEO query)"),
    ("1b",            "phase1b_curate_datasets.py",           "Dataset Curation (21 to 7 datasets)"),
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
    return [p[0] for p in PHASES]


def get_phase_by_id(phase_id):
    for pid, script, desc in PHASES:
        if pid == phase_id:
            return pid, script, desc
    return None


def format_duration(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    secs = seconds % 60
    if minutes < 60:
        return f"{minutes}m {secs:.0f}s"
    hours = int(minutes // 60)
    mins = minutes % 60
    return f"{hours}h {mins}m {secs:.0f}s"


def run_phase(phase_id, script, description, phase_num, total, log_dir=None):
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
            timeout=7200,
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
    return status, duration


def main():
    parser = argparse.ArgumentParser(
        description="Cross-Disease Transcriptomic Meta-Analysis Pipeline Runner",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="Run all pipeline phases")
    group.add_argument("--phase", nargs="+", metavar="ID", help="Run specific phase(s)")
    group.add_argument("--from-phase", metavar="ID", help="Run from specified phase")
    group.add_argument("--list", action="store_true", help="List all phases")

    parser.add_argument("--log", action="store_true", help="Save output to log files")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop on first failure")

    args = parser.parse_args()

    if args.list:
        for phase_id, script, desc in PHASES:
            print(f"  {phase_id:<15} {script:<42} {desc}")
        return

    if args.all:
        phases_to_run = PHASES[:]
    elif args.phase:
        phases_to_run = []
        for pid in args.phase:
            phase = get_phase_by_id(pid)
            if phase is None:
                print(f"Error: Unknown phase '{pid}'")
                sys.exit(1)
            phases_to_run.append(phase)
    elif args.from_phase:
        phase_ids = get_phase_ids()
        if args.from_phase not in phase_ids:
            print(f"Error: Unknown phase '{args.from_phase}'")
            sys.exit(1)
        start_idx = phase_ids.index(args.from_phase)
        phases_to_run = PHASES[start_idx:]

    log_dir = LOGS_DIR if args.log else None
    results = []
    total = len(phases_to_run)

    for i, (phase_id, script, desc) in enumerate(phases_to_run, 1):
        status, duration = run_phase(phase_id, script, desc, i, total, log_dir)
        results.append((phase_id, desc, status, duration))
        if status != 0 and args.stop_on_error:
            break

    failed = sum(1 for _, _, s, _ in results if s != 0)
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
