#!/usr/bin/env python3
"""
Phase 1: Systematic GEO Dataset Discovery for SpA Cross-Disease Meta-Analysis
Queries NCBI E-utilities to find relevant transcriptomic datasets.
"""

import requests
import json
import csv
import time
import os
from datetime import datetime

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUTPUT_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa/data/metadata"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SEARCH_QUERIES = {
    "AS": [
        '"ankylosing spondylitis"[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC" OR "whole blood")',
    ],
    "PsA": [
        '"psoriatic arthritis"[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC")',
    ],
    "IBD_arthritis": [
        '("inflammatory bowel disease" OR "Crohn")[Title] AND "Homo sapiens"[Organism] AND ("blood" OR "PBMC")',
    ],
}


def search_geo(query, db="gds", retmax=100):
    url = f"{BASE_URL}/esearch.fcgi"
    params = {"db": db, "term": query, "retmax": retmax, "retmode": "json"}
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        print(f"  Found {len(ids)} results")
        return ids
    except Exception as e:
        print(f"  ERROR: {e}")
        return []


def fetch_geo_details(id_list, db="gds"):
    if not id_list:
        return []
    results = []
    for i in range(0, len(id_list), 20):
        batch = id_list[i:i+20]
        url = f"{BASE_URL}/esummary.fcgi"
        params = {"db": db, "id": ",".join(batch), "retmode": "json"}
        try:
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            result_data = data.get("result", {})
            for uid in batch:
                if uid in result_data:
                    results.append(result_data[uid])
        except Exception as e:
            print(f"  ERROR: {e}")
        time.sleep(0.5)
    return results


def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    all_results = {}

    print("Phase 1: GEO Dataset Discovery")
    for disease, queries in SEARCH_QUERIES.items():
        disease_ids = set()
        for q in queries:
            ids = search_geo(q)
            disease_ids.update(ids)
            time.sleep(1)
        details = fetch_geo_details(list(disease_ids))
        all_results[disease] = details

    raw_path = os.path.join(OUTPUT_DIR, f"geo_search_raw_{timestamp}.json")
    with open(raw_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Results saved to: {raw_path}")


if __name__ == "__main__":
    main()
