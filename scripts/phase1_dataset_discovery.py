#!/usr/bin/env python3
"""
Phase 1: Systematic GEO Dataset Discovery for SpA Cross-Disease Meta-Analysis
Queries NCBI E-utilities to find relevant transcriptomic datasets.
"""

import requests
import xml.etree.ElementTree as ET
import json
import csv
import time
import os
from datetime import datetime

# Configuration
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUTPUT_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa/data/metadata"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Search queries for each disease category
SEARCH_QUERIES = {
    "AS": [
        '"ankylosing spondylitis"[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC" OR "whole blood" OR "mononuclear")',
        '"ankylosing spondylitis"[Description] AND "Homo sapiens"[Organism] AND "expression profiling"[DataSet Type]',
    ],
    "PsA": [
        '"psoriatic arthritis"[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC" OR "whole blood")',
        '"psoriatic arthritis"[Description] AND "Homo sapiens"[Organism] AND "expression profiling"[DataSet Type]',
    ],
    "SpA_general": [
        '"spondyloarthritis"[Title] AND "Homo sapiens"[Organism]',
        '"spondyloarthropathy"[Title] AND "Homo sapiens"[Organism]',
    ],
    "ReA": [
        '"reactive arthritis"[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC")',
    ],
    "IBD_arthritis": [
        '("inflammatory bowel disease" OR "Crohn" OR "ulcerative colitis")[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC" OR "whole blood") AND "expression profiling"[DataSet Type]',
    ],
    "Uveitis": [
        '"uveitis"[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC")',
        '"anterior uveitis"[Title] AND "Homo sapiens"[Organism]',
    ],
}


def search_geo(query, db="gds", retmax=100):
    """Search GEO DataSets using E-utilities."""
    url = f"{BASE_URL}/esearch.fcgi"
    params = {
        "db": db,
        "term": query,
        "retmax": retmax,
        "retmode": "json",
        "usehistory": "y",
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        result = data.get("esearchresult", {})
        ids = result.get("idlist", [])
        count = result.get("count", "0")
        print(f"  Found {count} results, retrieved {len(ids)} IDs")
        return ids
    except Exception as e:
        print(f"  ERROR in search: {e}")
        return []


def fetch_geo_details(id_list, db="gds"):
    """Fetch detailed info for GEO DataSet IDs."""
    if not id_list:
        return []
    
    results = []
    # Process in batches of 20
    for i in range(0, len(id_list), 20):
        batch = id_list[i:i+20]
        url = f"{BASE_URL}/esummary.fcgi"
        params = {
            "db": db,
            "id": ",".join(batch),
            "retmode": "json",
        }
        try:
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            result_data = data.get("result", {})
            
            for uid in batch:
                if uid in result_data:
                    entry = result_data[uid]
                    results.append(entry)
        except Exception as e:
            print(f"  ERROR fetching batch: {e}")
        
        time.sleep(0.5)  # Rate limiting
    
    return results


def parse_geo_entry(entry):
    """Extract relevant fields from a GEO entry."""
    accession = entry.get("accession", "")
    title = entry.get("title", "")
    summary = entry.get("summary", "")
    gpl = entry.get("gpl", "")
    gse = entry.get("gse", "")
    taxon = entry.get("taxon", "")
    entrytype = entry.get("entrytype", "")
    n_samples = entry.get("n_samples", "")
    
    # Get sample counts
    samples = entry.get("samples", [])
    if isinstance(samples, list):
        n_samples = len(samples)
    
    # Try to determine if it's GSE
    if gse:
        accession = f"GSE{gse}" if not str(gse).startswith("GSE") else gse
    
    return {
        "uid": entry.get("uid", ""),
        "accession": accession,
        "title": title,
        "summary": summary[:500],  # Truncate
        "gpl": gpl,
        "taxon": taxon,
        "entrytype": entrytype,
        "n_samples": n_samples,
    }


def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    all_results = {}
    all_accessions = set()
    
    print("=" * 80)
    print("Phase 1: Systematic GEO Dataset Discovery for SpA Meta-Analysis")
    print(f"Timestamp: {timestamp}")
    print("=" * 80)
    
    for disease, queries in SEARCH_QUERIES.items():
        print(f"\n{'='*60}")
        print(f"Searching for: {disease}")
        print(f"{'='*60}")
        
        disease_ids = set()
        for q in queries:
            print(f"\nQuery: {q[:100]}...")
            ids = search_geo(q)
            disease_ids.update(ids)
            time.sleep(1)  # Rate limiting between queries
        
        print(f"\nTotal unique IDs for {disease}: {len(disease_ids)}")
        
        if disease_ids:
            details = fetch_geo_details(list(disease_ids))
            parsed = [parse_geo_entry(d) for d in details]
            all_results[disease] = parsed
            
            for p in parsed:
                if p["accession"]:
                    all_accessions.add(p["accession"])
            
            print(f"Parsed {len(parsed)} entries for {disease}")
        else:
            all_results[disease] = []
        
        time.sleep(1)
    
    # Save raw results
    raw_path = os.path.join(OUTPUT_DIR, f"geo_search_raw_{timestamp}.json")
    with open(raw_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nRaw results saved to: {raw_path}")
    
    # Create summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    for disease, entries in all_results.items():
        print(f"\n{disease}: {len(entries)} datasets found")
        for e in entries[:10]:
            print(f"  {e['accession']:15s} | {e['n_samples']:>4s} samples | {e['title'][:80]}")
    
    print(f"\nTotal unique accessions: {len(all_accessions)}")
    
    # Save summary CSV
    csv_path = os.path.join(OUTPUT_DIR, "geo_search_summary.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["disease_category", "accession", "title", "n_samples", "platform", "entrytype"])
        for disease, entries in all_results.items():
            for e in entries:
                writer.writerow([
                    disease, e["accession"], e["title"][:200], 
                    e["n_samples"], e["gpl"], e["entrytype"]
                ])
    print(f"Summary CSV saved to: {csv_path}")


if __name__ == "__main__":
    main()
