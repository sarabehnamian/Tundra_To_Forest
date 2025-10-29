#!/usr/bin/env python3
# 00_neotoma_download_extract_geo.py  (ALL POLLEN DATASETS, robust enumeration or numeric-ID range)
# Enumerate pollen datasets via /data/sites?datasettype=pollen (with fallbacks),
# OR scan a numeric ID range (for SLURM arrays),
# download bundles, extract ALL pollen (and optional spores) rows,
# add AGE metadata (including calendar-year columns), write neotoma_out/pollen_records.txt

import json, time, random, argparse, requests, pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Iterable

OUTDIR = Path("neotoma_out")

BASES = [
    "https://api.neotomadb.org/v2.0",
    "https://api.neotomadb.org/v1",
    "https://api-dev.neotomadb.org/v2.0",
]

S = requests.Session()
S.trust_env = False
S.headers.update({
    "User-Agent": "Neotoma pollen downloader-extractor (ages to CE)",
    "Accept": "application/json",
    "Connection": "close",
})

# ------------------ HTTP helpers ------------------
def backoff_get(url: str, params: Optional[Dict]=None, tries: int=7, timeout: int=60, debug: bool=False):
    last = None
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if debug:
                print(f"GET {r.url} -> {r.status_code}")
            if r.status_code >= 500:
                last = RuntimeError(f"{r.status_code} @ {r.url}")
                time.sleep(min(2**i, 20) + random.random()); continue
            r.raise_for_status()
            ct = r.headers.get("Content-Type","")
            if "json" in ct.lower():
                return r.json()
            try:
                return json.loads(r.text)
            except Exception:
                return {"data": None}
        except requests.RequestException as e:
            last = e
            time.sleep(min(2**i, 20) + random.random())
    if last:
        raise last
    raise RuntimeError("unreachable")

def get_json_any(path: str, params: Optional[Dict]=None, debug: bool=False):
    err = None
    for base in BASES:
        try:
            return backoff_get(f"{base}{path}", params=params, debug=debug)
        except Exception as e:
            err = e
            if debug:
                print(f"[WARN] Base failed: {base}{path} :: {e}")
            continue
    raise RuntimeError(f"All Neotoma bases failed for {path}: {err}")

# ------------------ utilities ------------------
def ensure_list(x: Any) -> List[Any]:
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]

def parse_geojson_string(s: Any) -> Tuple[Optional[float], Optional[float]]:
    try:
        gj = json.loads(s) if isinstance(s, str) else s
        if isinstance(gj, dict) and "coordinates" in gj:
            lon, lat = gj["coordinates"][:2]
            return float(lat), float(lon)
    except Exception:
        pass
    return None, None

def get_default_chronology_id(collectionunit: Dict) -> Optional[int]:
    dc = collectionunit.get("defaultchronology")
    if isinstance(dc, int):
        return dc
    if isinstance(dc, dict):
        cid = dc.get("chronologyid")
        try:
            return int(cid) if cid is not None else None
        except Exception:
            return None
    return None

def bp_to_ce(val: Optional[float]) -> Optional[float]:
    """Convert years BP (Before Present, 1950) to calendar year CE (negative == BCE)."""
    try:
        if val is None:
            return None
        return 1950.0 - float(val)
    except Exception:
        return None

def pick_age_record(sample: Dict, default_chron_id: Optional[int]) -> Tuple[Optional[Dict], str]:
    ages = (sample.get("ages") or [])
    if default_chron_id is not None:
        for a in ages:
            if isinstance(a, dict) and a.get("age") is not None and a.get("chronologyid") == default_chron_id:
                return a, "default_chronology"
    for a in ages:
        if isinstance(a, dict) and a.get("age") is not None:
            return a, "first_nonnull_age"
    if any(sample.get(k) is not None for k in ("age","ageolder","ageyounger")):
        return {
            "age": sample.get("age"),
            "ageolder": sample.get("ageolder"),
            "ageyounger": sample.get("ageyounger"),
            "agetype": sample.get("agetype"),
            "chronologyid": sample.get("chronologyid"),
            "chronologyname": sample.get("chronologyname"),
        }, "sample_fields"
    return None, "none"

# ------------------ dataset enumeration via SITES ------------------
def _collect_ids_from_sites_payload(payload: Dict, kw_lower: str, debug: bool=False) -> List[int]:
    """Walk the /data/sites payload to extract dataset IDs of type containing 'pollen'."""
    data = payload.get("data", payload)
    sites: Iterable = data if isinstance(data, list) else [data]
    out: List[int] = []
    n_sites = 0
    for site in sites:
        if not isinstance(site, dict):
            continue
        n_sites += 1
        collunits = (site.get("collunits") or site.get("collectionunits") or site.get("collectionunit") or [])
        for cu in ensure_list(collunits):
            dsets = (cu.get("datasets") or cu.get("dataset") or [])
            for ds in ensure_list(dsets):
                if not isinstance(ds, dict):
                    continue
                dtype = (ds.get("datasettype") or ds.get("DatasetType") or "")
                if kw_lower in str(dtype).lower():
                    dsid = ds.get("datasetid") or ds.get("DatasetID") or ds.get("datasetId")
                    try:
                        if dsid is not None:
                            out.append(int(dsid))
                    except Exception:
                        pass
    if debug:
        print(f"[DEBUG] Sites seen: {n_sites}, dataset IDs collected this pass: {len(out)}")
    return out

def fetch_pollen_dataset_ids(limit: int = 200, keyword: str = "pollen", try_all_data: bool = True, debug: bool=False) -> List[int]:
    """
    Enumerate ALL dataset IDs by querying /data/sites?datasettype=pollen (server-side filter),
    walking nested collunits/datasets, and collecting datasetids whose datasettype contains `keyword`.
    Tries all_data=true first; falls back to pagination.
    """
    kw = (keyword or "").lower()
    ids: List[int] = []

    # First attempt: all_data=true
    if try_all_data:
        try:
            j = get_json_any("/data/sites", params={"datasettype": "pollen", "all_data": "true"}, debug=debug)
            ids = _collect_ids_from_sites_payload(j, kw, debug=debug)
            if ids:
                return sorted(set(ids))
        except Exception as e:
            if debug:
                print(f"[DEBUG] all_data=true failed: {e}")

    # Paginated fallback
    offset = 0
    empty_hits = 0
    while True:
        params = {"datasettype": "pollen", "limit": limit, "offset": offset}
        try:
            j = get_json_any("/data/sites", params=params, debug=debug)
        except Exception as e:
            if debug:
                print(f"[DEBUG] pagination fetch failed at offset={offset}: {e}")
            break
        new_ids = _collect_ids_from_sites_payload(j, kw, debug=debug)
        if not new_ids:
            empty_hits += 1
        else:
            ids.extend(new_ids)
            empty_hits = 0
        offset += limit
        if empty_hits >= 2:
            break
        time.sleep(0.05)

    return sorted(set(ids))

# ------------------ extraction ------------------
def extract_bundle_rows(bundle: Dict, dataset_id: int, include_spores: bool) -> List[Dict]:
    rows: List[Dict] = []

    site = bundle.get("site", {}) or {}
    sitename = site.get("sitename")
    altitude = site.get("altitude")
    lat, lon = parse_geojson_string(site.get("geography"))

    cu = site.get("collectionunit", {}) or {}
    ds = cu.get("dataset", {}) or site.get("dataset", {}) or {}
    dset_type = ds.get("datasettype")
    coll_date = cu.get("colldate") or ds.get("collectiondate")

    default_chron_id = get_default_chronology_id(cu)

    for samp in ds.get("samples", []) or []:
        age_rec, age_src = pick_age_record(samp, default_chron_id)
        age_bp = age_old = age_young = age_type = chrono_name = None
        chrono_id = None
        if age_rec:
            age_bp    = age_rec.get("age")
            age_old   = age_rec.get("ageolder")
            age_young = age_rec.get("ageyounger")
            age_type  = age_rec.get("agetype")
            chrono_name = age_rec.get("chronologyname")
            chrono_id   = age_rec.get("chronologyid")

        sample_id = samp.get("sampleid")
        depth_cm = samp.get("depth")

        taxon_items = samp.get("datum", []) or samp.get("data", []) or []
        for d in taxon_items:
            element = (d.get("element") or d.get("elementtype") or "").lower()
            if ("pollen" in element) or (include_spores and "spore" in element):
                rows.append({
                    "dataset_id": dataset_id,
                    "site_name": sitename,
                    "latitude": lat,
                    "longitude": lon,
                    "altitude": altitude,
                    "dataset_type": dset_type,
                    "collection_date": coll_date,
                    "sample_id": sample_id,
                    "depth_cm": depth_cm,

                    # AGE FIELDS (BP & metadata)
                    "age_BP": age_bp,
                    "date_BP": age_bp,  # alias
                    "age_older": age_old,
                    "age_younger": age_young,
                    "age_type": age_type,
                    "chronology_name": chrono_name,
                    "chronology_id": chrono_id,
                    "age_source": age_src,

                    # TAXON DATA
                    "taxon": d.get("variablename") or d.get("taxonname"),
                    "count": d.get("value"),
                    "units": d.get("units"),
                    "element": element,
                })
    return rows

# ------------------ I/O helpers ------------------
def read_ids(ids_path: Optional[str], id_list: Optional[str]) -> List[int]:
    ids: List[int] = []
    if ids_path:
        p = Path(ids_path)
        if p.exists():
            for line in p.read_text(encoding="utf-8").splitlines():
                t = line.strip()
                if t.isdigit():
                    ids.append(int(t))
    if id_list:
        for t in id_list.split(","):
            t = t.strip()
            if t.isdigit():
                ids.append(int(t))
    return sorted(set(ids))

# ------------------ main ------------------
def main():
    ap = argparse.ArgumentParser(description="Download ALL pollen datasets from Neotoma and extract ALL pollen rows with ages (incl. CE columns).")
    ap.add_argument("--ids", help="Path to a file with dataset IDs (one per line).")
    ap.add_argument("--id-list", help="Comma-separated IDs to add (e.g., 2637,3656).")
    ap.add_argument("--id-start", type=int, help="Start dataset ID (inclusive) for numeric range scan (for SLURM arrays).")
    ap.add_argument("--id-end",   type=int, help="End dataset ID (inclusive) for numeric range scan (for SLURM arrays).")
    ap.add_argument("--max-ids", type=int, default=None, help="Cap total IDs after enumeration (testing).")
    ap.add_argument("--outdir", default=str(OUTDIR), help="Output folder (default: neotoma_out)")
    ap.add_argument("--no-spores", action="store_true", help="Exclude spores (pollen only).")
    ap.add_argument("--throttle", type=float, default=0.4, help="Sleep seconds between downloads.")
    ap.add_argument("--debug", action="store_true", help="Verbose fetch logs for troubleshooting.")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    include_spores = not args.no_spores

    # --- Choose IDs: either a numeric range (for SLURM arrays) OR enumerate all pollen datasets
    user_ids = read_ids(args.ids, args.id_list)
    ids: List[int] = []

    if (args.id_start is not None) or (args.id_end is not None):
        if (args.id_start is None) or (args.id_end is None):
            raise SystemExit("ERROR: Use --id-start and --id-end together.")
        if args.id_start > args.id_end:
            raise SystemExit("ERROR: --id-start must be <= --id-end.")
        ids = list(range(args.id_start, args.id_end + 1))
        print(f"Using explicit numeric ID range: {args.id_start}-{args.id_end}")
    else:
        print("Enumerating pollen dataset IDs via /data/sites?datasettype=pollen ...")
        enum_ids = fetch_pollen_dataset_ids(limit=200, keyword="pollen", try_all_data=True, debug=args.debug)
        if args.max_ids is not None:
            enum_ids = enum_ids[:args.max_ids]
        print(f"Found {len(enum_ids)} pollen dataset IDs via sites enumeration.")
        ids = enum_ids

    if user_ids:
        print(f"Adding {len(user_ids)} user-specified IDs.")
    ids = sorted(set(ids).union(user_ids))

    print(f"Total unique dataset IDs to process: {len(ids)}")
    if not ids:
        print("No rows extracted.")
        return

    all_rows: List[Dict] = []

    for i, did in enumerate(ids, 1):
        raw_path = outdir / f"raw_{did}.json"
        j = None
        src = "downloaded"
        if raw_path.exists():
            try:
                j = json.loads(raw_path.read_text(encoding="utf-8"))
                src = "cached"
            except Exception:
                src = "cache_read_error"

        if j is None:
            try:
                j = get_json_any(f"/data/downloads/{did}", debug=args.debug)
            except Exception as e:
                print(f"[{i}/{len(ids)}] {did} ERROR: {e}")
                time.sleep(args.throttle)
                continue
            try:
                raw_path.write_text(json.dumps(j, indent=2), encoding="utf-8")
            except Exception:
                pass

        root = j.get("data", j)
        bundles = root if isinstance(root, list) else [root]
        rows_here = []
        for b in bundles:
            rows_here.extend(extract_bundle_rows(b, did, include_spores))
        print(f"[{i}/{len(ids)}] {did} ({src}) -> +{len(rows_here)} rows")
        all_rows.extend(rows_here)
        time.sleep(args.throttle)

    if not all_rows:
        print("No rows extracted.")
        return

    # === Calendar-year columns (CE; negative == BCE) ===
    for r in all_rows:
        r["year_CE"]         = bp_to_ce(r.get("age_BP"))
        r["year_older_CE"]   = bp_to_ce(r.get("age_older"))
        r["year_younger_CE"] = bp_to_ce(r.get("age_younger"))
        # Midpoint CE if a range exists; otherwise fall back to single-point age
        if (r["year_older_CE"] is not None) and (r["year_younger_CE"] is not None):
            r["year_midpoint_CE"] = (r["year_older_CE"] + r["year_younger_CE"]) / 2.0
        else:
            r["year_midpoint_CE"] = r["year_CE"]

    out_csv = outdir / "pollen_records.txt"  # CSV-in-.txt
    pd.DataFrame(all_rows).to_csv(out_csv, index=False)
    print(f"\nWrote {len(all_rows)} rows -> {out_csv}")

if __name__ == "__main__":
    main()
