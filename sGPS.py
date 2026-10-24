#!/usr/bin/env python3
# sGPS.py — GIGGLE-like index/search (BGZF+tabix, sorting, bp/bins Fisher, two types of scores,
# chr-normalization, bin-cache + precomputed cache, two-sided Fisher option, FDR (BH))
# CLI:
#   INDEX:
#     python sGPS.py index <db_dir> --genome genome.txt [-t 4]
#       [--precompute-bins] [--bin-size 4096]
#   SEARCH:
#     python sGPS.py search <query.bed[.gz]> <db_dir> --genome genome.txt \
#       [--stat {bins,bp}] [--score {giggle,odds}] [--bin-size 4096] \
#       [--normalize-chr] [--cache-bins] [--fisher {one,two}] [--fdr] \
#       [-t 4] [-o out.tsv]

import argparse, os, sys, json, shutil, math, gzip, tempfile, hashlib, subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.stats import fisher_exact
import pysam

MANIFEST = "sGPS_manifest.json"
META_DIR = ".sGPS"
CACHE_DIR = "cache"
META_VERSION = 5  # bump when the manifest or cache format changes

def die(msg, code=1):
    print(f"[sGPS] ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

# ---------- genome / names ----------

def read_genome(genome_file):
    sizes = {}
    total = 0
    with open(genome_file) as fh:
        for line in fh:
            if not line.strip(): continue
            chrom, length = line.strip().split()[:2]
            l = int(length)
            sizes[chrom] = l
            total += l
    if not sizes:
        die("Empty genome file.")
    return sizes, total

def genome_has_chr_prefix(chrom_sizes):
    keys = list(chrom_sizes.keys())
    if not keys: return False
    pref = sum(1 for k in keys if k.startswith("chr"))
    return pref >= (len(keys) / 2.0)

def add_chr(name):   return name if name.startswith("chr") else "chr"+name
def strip_chr(name): return name[3:] if name.startswith("chr") else name

def build_name_mapper(chrom_sizes, normalize_chr):
    if not normalize_chr:
        return (lambda x: x)
    want_chr = genome_has_chr_prefix(chrom_sizes)
    return add_chr if want_chr else strip_chr

# ---------- manifest ----------

def ensure_meta_dir(db_dir):
    meta_root = Path(db_dir)/META_DIR
    meta_root.mkdir(exist_ok=True, parents=True)
    (meta_root/ CACHE_DIR).mkdir(exist_ok=True, parents=True)
    return meta_root

def load_manifest(db_dir):
    mf = Path(db_dir)/META_DIR/MANIFEST
    if not mf.exists(): return None
    with open(mf) as fh:
        return json.load(fh)

def save_manifest(db_dir, manifest):
    mf = Path(db_dir)/META_DIR/MANIFEST
    with open(mf, "w") as fh:
        json.dump(manifest, fh, indent=2)

def list_bed_files(db_dir):
    exts = (".bed", ".bed.gz")
    for root, _, files in os.walk(db_dir):
        for f in files:
            if f.lower().endswith(exts):
                yield Path(root)/f

# ---------- sorting / bgzf / tabix ----------

def _is_tabix_openable(path):
    try:
        pysam.TabixFile(str(path)).close()
        return True
    except Exception:
        return False

def _python_sort_bed(src_path, out_path):
    header = []
    rows = []
    with open(src_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#","track","browser")):
                header.append(line); continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                try:
                    rows.append((parts[0], int(parts[1]), line))
                except:
                    header.append(line)
            else:
                header.append(line)
    with open(out_path, "w") as out:
        for h in header: out.write(h)
        for _, __, line in sorted(rows, key=lambda x: (x[0], x[1])):
            out.write(line)

def _sort_bed_to_tmp(src_path):
    tmp_fd, tmp_name = tempfile.mkstemp(suffix=".bed")
    os.close(tmp_fd)
    tmp = Path(tmp_name)
    try:
        env = dict(os.environ)
        env["LC_ALL"] = "C"
        with open(tmp, "w") as out:
            subprocess.run(["sort", "-k1,1", "-k2,2n", str(src_path)],
                           check=True, stdout=out, env=env)
        return tmp
    except Exception:
        _python_sort_bed(src_path, tmp)
        return tmp

def _compress_and_index(plain_path, preset="bed"):
    gz_path = Path(str(plain_path)+".gz")
    pysam.tabix_compress(str(plain_path), str(gz_path), force=True)
    pysam.tabix_index(str(gz_path), preset=preset, force=True)
    if not _is_tabix_openable(gz_path):
        tmp_sorted = _sort_bed_to_tmp(plain_path)
        pysam.tabix_compress(str(tmp_sorted), str(gz_path), force=True)
        pysam.tabix_index(str(gz_path), preset=preset, force=True)
        os.remove(tmp_sorted)
    return gz_path, Path(str(gz_path)+".tbi")

def bgzip_tabix_if_needed(bed_path, preset="bed"):
    bed_path = Path(bed_path)
    if bed_path.suffix != ".gz":
        return _compress_and_index(bed_path, preset=preset)
    gz_path = bed_path
    tbi_path = Path(str(gz_path)+".tbi")
    try:
        if not tbi_path.exists():
            pysam.tabix_index(str(gz_path), preset=preset, force=True)
        if _is_tabix_openable(gz_path):
            return gz_path, tbi_path
    except Exception:
        pass
    bgzf_path = Path(str(bed_path)+".bgzf.gz")
    with gzip.open(bed_path, "rb") as fin, tempfile.NamedTemporaryFile(delete=False) as tmp_plain:
        shutil.copyfileobj(fin, tmp_plain)
        tmp_plain_path = Path(tmp_plain.name)
    tmp_sorted = _sort_bed_to_tmp(tmp_plain_path)
    pysam.tabix_compress(str(tmp_sorted), str(bgzf_path), force=True)
    pysam.tabix_index(str(bgzf_path), preset=preset, force=True)
    os.remove(tmp_plain_path); os.remove(tmp_sorted)
    return bgzf_path, Path(str(bgzf_path)+".tbi")

# ---------- intervals / coverage ----------

def merge_intervals(lst):
    if not lst: return []
    lst.sort()
    out=[]; s,e=lst[0]
    for x,y in lst[1:]:
        if x <= e: e = max(e,y)
        else: out.append((s,e)); s,e = x,y
    out.append((s,e))
    return out

def merge_coverage_bp(intervals):
    if not intervals: return 0
    merged = merge_intervals(intervals[:])
    return sum(e - s for s, e in merged)

def compute_file_coverage_bp_and_count(tabix_path, genome_sizes):
    tb = pysam.TabixFile(str(tabix_path))
    total_bp = 0
    interval_count = 0
    for chrom in genome_sizes.keys():
        try:
            ivals = []
            for rec in tb.fetch(chrom):
                parts = rec.split("\t")
                if len(parts) < 3: continue
                start = int(parts[1]); end = int(parts[2])
                if end > start:
                    ivals.append((start,end))
                    interval_count += 1
            total_bp += merge_coverage_bp(ivals)
        except ValueError:
            continue
    tb.close()
    return total_bp, interval_count

# ---------- bins & cache ----------

def build_bins_for_genome(genome_sizes, bin_size):
    chrom_bins_len = {c: ((genome_sizes[c] + bin_size - 1)//bin_size) for c in genome_sizes}
    total_bins = sum(chrom_bins_len.values())
    return chrom_bins_len, total_bins

def mark_bins_from_intervals(bin_arr, intervals, bin_size):
    n = bin_arr.shape[0]
    for s,e in intervals:
        if e <= s: continue
        b0 = max(0, min(n-1, s//bin_size))
        b1 = max(0, min(n-1, (e-1)//bin_size))
        if b1 >= b0:
            bin_arr[b0:b1+1] = True

def iter_db_intervals(tb, chrom, start=None, end=None):
    try:
        it = tb.fetch(chrom) if start is None else tb.fetch(chrom, start, end)
        for rec in it:
            parts = rec.split("\t")
            if len(parts) < 3: continue
            rs = int(parts[1]); re = int(parts[2])
            if re > rs: yield (rs,re)
    except ValueError:
        return

def _fisher_table_counts_from_bins(q_bins, d_bins, total_bins):
    both  = np.logical_and(q_bins,  d_bins).sum()
    q_only= np.logical_and(q_bins, ~d_bins).sum()
    d_only= np.logical_and(~q_bins, d_bins).sum()
    neither = int(total_bins) - int(both + q_only + d_only)
    return int(both), int(q_only), int(d_only), int(neither)

def fisher_from_bins(q_bins, d_bins, total_bins, fisher_mode="one"):
    both, q_only, d_only, neither = _fisher_table_counts_from_bins(q_bins, d_bins, total_bins)
    if (both+q_only)==0 or (d_only+neither)==0 or (both+d_only)==0 or (q_only+neither)==0:
        odds, pval = 0.0, 1.0
    else:
        alt = 'greater' if fisher_mode == "one" else 'two-sided'
        odds, pval = fisher_exact([[both, q_only],[d_only, neither]], alternative=alt)
        if np.isinf(odds):
            odds = float((both*(neither+1))/((q_only+1)*max(d_only,1)))
    return both, q_only, d_only, neither, float(odds), float(pval)

def fisher_from_bp(a_overlap_bp, q_bp, d_bp, genome_bp, fisher_mode="one"):
    a = int(a_overlap_bp)
    b = max(int(q_bp - a), 0)
    c = max(int(d_bp - a), 0)
    d = max(int(genome_bp - q_bp - d_bp + a), 0)
    if (a+b)==0 or (c+d)==0 or (a+c)==0 or (b+d)==0:
        odds, pval = 0.0, 1.0
    else:
        alt = 'greater' if fisher_mode == "one" else 'two-sided'
        odds, pval = fisher_exact([[a,b],[c,d]], alternative=alt)
        if np.isinf(odds):
            odds = float((a*(d+1))/((b+1)*max(c,1)))
    return float(odds), float(pval)

def combo_score(odds, pval, mode="giggle"):
    pval = max(min(float(pval), 1.0), 1e-300)
    if odds <= 0:
        return 0.0
    if mode == "giggle":
        return (-math.log10(pval)) * (math.log2(odds))
    else:
        return (-math.log10(pval)) * odds

def overlap_bp_of_lists(a_list, b_list):
    i=j=0; total=0
    a=sorted(a_list); b=sorted(b_list)
    while i < len(a) and j < len(b):
        s1,e1=a[i]; s2,e2=b[j]
        s=max(s1,s2); e=min(e1,e2)
        if e > s: total += (e - s)
        if e1 < e2: i+=1
        else: j+=1
    return total

def _cache_key_for_file(fmeta, bin_size, genome_sizes):
    h = hashlib.sha1()
    payload = json.dumps({
        "path": fmeta["bgz"],
        "mtime": fmeta.get("mtime"),
        "size": fmeta.get("size"),
        "bin_size": bin_size,
        "genome": genome_sizes,
        "ver": META_VERSION
    }, sort_keys=True).encode()
    h.update(payload)
    return h.hexdigest()

def _cache_path(db_dir, key):
    return Path(db_dir)/META_DIR/CACHE_DIR/(key+".npy")

def _load_bins_cache(db_dir, key):
    path = _cache_path(db_dir, key)
    if path.exists():
        try:
            arr = np.load(path, allow_pickle=False, mmap_mode='r')
            return arr.astype(bool, copy=False)
        except Exception:
            return None
    return None

def _save_bins_cache(db_dir, key, arr):
    path = _cache_path(db_dir, key)
    try:
        np.save(path, arr.astype(np.bool_), allow_pickle=False)
    except Exception:
        pass

# ---------- indexing ----------

def index_db(db_dir, genome_file, threads=4, precompute_bins=False, pre_bin_size=4096):
    genome_sizes, genome_bp = read_genome(genome_file)
    ensure_meta_dir(db_dir)

    files = list(list_bed_files(db_dir))
    if not files:
        die("No BED/BED.GZ files found under db_dir.")

    print(f"[sGPS] Indexing {len(files)} file(s) ...")

    manifest = {
        "version": META_VERSION,
        "genome_file": str(Path(genome_file).resolve()),
        "genome_total_bp": genome_bp,
        "chrom_sizes": genome_sizes,
        "files": [],
        "precomputed_bins": {"enabled": bool(precompute_bins),
                             "bin_size": int(pre_bin_size) if precompute_bins else None}
    }

    def work(path):
        gz_path, _ = bgzip_tabix_if_needed(path)
        cov_bp, n_int = compute_file_coverage_bp_and_count(gz_path, genome_sizes)
        return {
            "source": str(Path(path).resolve()),
            "bgz": str(Path(gz_path).resolve()),
            "coverage_bp": int(cov_bp),
            "intervals": int(n_int),
            "mtime": os.path.getmtime(path),
            "size": os.path.getsize(path),
        }

    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {ex.submit(work, p): p for p in files}
        for fut in as_completed(futs):
            info = fut.result()
            manifest["files"].append(info)
            print(f"[sGPS] Indexed: {Path(info['source']).name}  bp={info['coverage_bp']}  n={info['intervals']}")

    if precompute_bins:
        print(f"[sGPS] Precomputing bin masks (bin_size={pre_bin_size}) ...")
        chrom_bins_len, total_bins = build_bins_for_genome(genome_sizes, pre_bin_size)

        def build_and_save(fmeta):
            key = _cache_key_for_file(fmeta, pre_bin_size, genome_sizes)
            if _load_bins_cache(db_dir, key) is not None:
                return
            tb = pysam.TabixFile(str(fmeta["bgz"]))
            d_bins = np.zeros(total_bins, dtype=bool)
            offset = 0
            for chrom, n in chrom_bins_len.items():
                tmp=[]
                for rs,re in iter_db_intervals(tb, chrom):
                    tmp.append((rs,re))
                if tmp:
                    arr = np.zeros(n, dtype=bool)
                    mark_bins_from_intervals(arr, tmp, pre_bin_size)
                    d_bins[offset:offset+n] = arr
                offset += n
            tb.close()
            _save_bins_cache(db_dir, key, d_bins)

        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(build_and_save, fmeta) for fmeta in manifest["files"]]
            for _ in as_completed(futs):
                pass
        print("[sGPS] Precompute finished.")

    save_manifest(db_dir, manifest)
    print(f"[sGPS] Done. Manifest: {Path(db_dir)/META_DIR/MANIFEST}")

# ---------- query load ----------

def load_query(query_bed, name_map=lambda x:x):
    q = {}
    raw_count = 0
    opener = gzip.open if str(query_bed).endswith(".gz") else open
    with opener(query_bed, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#","track","browser")): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            chrom = name_map(parts[0])
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue
            if e <= s: continue
            q.setdefault(chrom, []).append((s,e))
            raw_count += 1
    q_merged_bp = sum(merge_coverage_bp(lst[:]) for lst in q.values())
    return q, q_merged_bp, raw_count

# ---------- FDR (BH) ----------

def bh_fdr(pvals):
    """Benjamini–Hochberg q-values. pvals: list[float]. Returns: list of q-values in the same order."""
    m = len(pvals)
    if m == 0:
        return []
    # Rank p-values
    order = np.argsort(pvals)
    ranks = np.empty(m, dtype=int)
    ranks[order] = np.arange(1, m+1)
    # BH formula: q_i = p_i * m / rank_i
    p = np.asarray(pvals, dtype=float)
    q = p * m / ranks
    # Enforce monotonicity from the end: q_i = min(q_i, q_{i+1})
    q_sorted = q[order]
    for i in range(m-2, -1, -1):
        if q_sorted[i] > q_sorted[i+1]:
            q_sorted[i] = q_sorted[i+1]
    q[order] = q_sorted
    # clamp [0,1]
    q = np.clip(q, 0.0, 1.0)
    return q.tolist()

# ---------- search ----------

def search(query_bed, db_dir, genome_file, out_tsv, threads=4, stat="bins",
           score_mode="giggle", bin_size=4096, normalize_chr=False, cache_bins=False,
           fisher_mode="two", do_fdr=False):
    manifest = load_manifest(db_dir)
    if manifest is None:
        die("No manifest found. Run 'index' first.")

    man_genome_file = manifest.get("genome_file")
    man_chrom_sizes = manifest.get("chrom_sizes", {})
    _, man_total_bp = read_genome(man_genome_file) if man_genome_file else (None, None)
    gsizes, genome_bp = read_genome(genome_file)

    if man_total_bp and man_total_bp != genome_bp:
        print("[sGPS] WARNING: genome total bp differs from manifest; using search-time genome.", file=sys.stderr)
    missing = set(gsizes.keys()) - set(man_chrom_sizes.keys())
    if missing:
        print(f"[sGPS] WARNING: chromosomes in --genome missing from manifest genome: {sorted(list(missing))}", file=sys.stderr)

    name_map = build_name_mapper(gsizes, normalize_chr)

    q_intervals, q_bp, q_raw_count = load_query(query_bed, name_map=name_map)
    if stat == "bp" and q_bp == 0:
        die("Query has zero merged coverage bp.")
    if stat == "bins" and q_raw_count == 0:
        die("Query has zero intervals.")

    files = manifest["files"]
    if not files: die("No indexed files in manifest.")

    # informational warnings
    query_chroms = set(q_intervals.keys())
    genome_chroms = set(gsizes.keys())
    missing_in_genome = sorted(query_chroms - genome_chroms)
    missing_in_query  = sorted(genome_chroms - query_chroms)
    if missing_in_genome:
        print(f"[sGPS] WARNING Query contains chromosomes not in genome file: {missing_in_genome}", file=sys.stderr)
    if missing_in_query:
        print(f"[sGPS] WARNING Genome contains chromosomes not present in query: {missing_in_query}", file=sys.stderr)

    for fmeta in files:
        try:
            tb = pysam.TabixFile(str(fmeta["bgz"]))
            db_chroms = set(tb.contigs)
            tb.close()
            missing_in_db = sorted(genome_chroms - db_chroms)
            extra_in_db   = sorted(db_chroms - genome_chroms)
            if missing_in_db or extra_in_db:
                print(
                    f"[sGPS] WARNING Chromosome mismatch in {os.path.basename(fmeta['source'])}: "
                    f"missing_in_db={missing_in_db}, extra_in_db={extra_in_db}",
                    file=sys.stderr
                )
        except Exception as e:
            print(f"[sGPS] WARNING Could not read contigs from {fmeta['source']}: {e}", file=sys.stderr)

    # bins preparation
    if stat == "bins":
        chrom_bins_len, total_bins = build_bins_for_genome(gsizes, bin_size)
        q_bins_concat=[]; chrom_offsets={}; offset=0
        for chrom in gsizes:
            n = chrom_bins_len[chrom]
            chrom_offsets[chrom]=offset
            arr = np.zeros(n, dtype=bool)
            if chrom in q_intervals:
                mark_bins_from_intervals(arr, q_intervals[chrom], bin_size)
            q_bins_concat.append(arr)
            offset += n
        q_bins = np.concatenate(q_bins_concat)
    else:
        chrom_bins_len = total_bins = q_bins = chrom_offsets = None

    # --- compute per-file statistics
    tmp_results=[]  # collecting here, before FDR calculation
    def work(fmeta):
        path = fmeta["bgz"]
        d_bp = int(fmeta["coverage_bp"])
        d_n  = int(fmeta["intervals"])
        tb = pysam.TabixFile(str(path))

        if stat == "bp":
            if d_bp == 0:
                tb.close()
                odds, pval, combo = 0.0, 1.0, 0.0
                return (fmeta["source"], 0, q_bp, d_bp, odds, pval, combo)
            total_overlap_bp = 0
            for chrom, qlist in q_intervals.items():
                if chrom not in gsizes:
                    continue
                q_merged = merge_intervals(qlist[:])
                tmp=[]
                for (qs,qe) in q_merged:
                    for rs,re in iter_db_intervals(tb, chrom, qs, qe):
                        tmp.append((rs,re))
                if not tmp: continue
                db_merged = merge_intervals(tmp)
                total_overlap_bp += overlap_bp_of_lists(q_merged, db_merged)
            odds, pval = fisher_from_bp(total_overlap_bp, q_bp, d_bp, genome_bp, fisher_mode=fisher_mode)
            tb.close()
            combo = combo_score(odds, pval, mode=score_mode)
            return (fmeta["source"], int(total_overlap_bp), int(q_bp), int(d_bp),
                    float(odds), float(pval), float(combo))

        else:
            # bins mode – use cache if available
            d_bins = None
            key = _cache_key_for_file(fmeta, bin_size, gsizes)
            if cache_bins:
                d_bins = _load_bins_cache(db_dir, key)
            if d_bins is None:
                d_bins = np.zeros(total_bins, dtype=bool)
                for chrom in gsizes.keys():
                    n = chrom_bins_len[chrom]
                    offs = chrom_offsets[chrom]
                    tmp=[]
                    for rs,re in iter_db_intervals(tb, chrom):
                        tmp.append((rs,re))
                    if tmp:
                        arr = np.zeros(n, dtype=bool)
                        mark_bins_from_intervals(arr, tmp, bin_size)
                        d_bins[offs:offs+n] = arr
                if cache_bins:
                    _save_bins_cache(db_dir, key, d_bins)
            tb.close()
            both, q_only, d_only, neither, odds, pval = fisher_from_bins(q_bins, d_bins, total_bins, fisher_mode=fisher_mode)
            combo = combo_score(odds, pval, mode=score_mode)
            return (fmeta["source"], int(both), int(q_raw_count), int(d_n),
                    float(odds), float(pval), float(combo))

    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [ex.submit(work, f) for f in files]
        for fut in as_completed(futs):
            tmp_results.append(fut.result())

    # FDR calculation (optional)
    if do_fdr:
        pvals = [row[-2] for row in tmp_results]  # Pvalue column in the rows
        qvals = bh_fdr(pvals)
        tmp_results = [tuple(list(row) + [float(qvals[i])]) for i, row in enumerate(tmp_results)]

    # sorting by combo score (q does not influence the order)
    results = sorted(tmp_results, key=lambda x: x[-2] if do_fdr else x[-1], reverse=True) if False else \
              sorted(tmp_results, key=lambda x: x[-2] if (do_fdr and isinstance(x[-1], float)) else x[-1], reverse=True)

    # header
    header_bp   = ["DB_File","OverlapBp","QueryBp","DBBp","OddsRatio","Pvalue","ComboScore"]
    header_bins = ["DB_File","OverlapBins","QueryIntervals","DBIntervals","OddsRatio","Pvalue","ComboScore"]
    header = header_bp if stat=="bp" else header_bins
    if do_fdr:
        header = header + ["Qvalue"]

    outfh = sys.stdout if out_tsv in ("-", "", None) else open(out_tsv, "w")
    print("\t".join(header), file=outfh)
    for row in results:
        print("\t".join([str(x) if x is not None else "" for x in row]), file=outfh)
    if outfh is not sys.stdout: outfh.close()

# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(prog="sGPS", description="GIGGLE-like genomic interval index & search (bins/BGZF/Fisher)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    api = sub.add_parser("index", help="Index all BED/BED.GZ files under a directory (BGZF + tabix + coverage).")
    api.add_argument("db_dir")
    api.add_argument("--genome", required=True, help="Genome sizes file (chrom<TAB>length), same build as the query and DB.")
    api.add_argument("-t","--threads", type=int, default=4)
    api.add_argument("--precompute-bins", action="store_true",
                     help="Precompute and save per-file bin masks to accelerate searches.")
    api.add_argument("--bin-size", type=int, default=4096,
                     help="Bin size in bp for precomputed bin masks. Use the same value for subsequent searches.")

    aps = sub.add_parser("search", help="Search query BED against indexed database.")
    aps.add_argument("query_bed")
    aps.add_argument("db_dir")
    aps.add_argument("--genome", required=True, help="Genome sizes file (chrom<TAB>length), same build as the query and DB.")
    aps.add_argument("-o","--out", default="-", help="Output TSV (default stdout).")
    aps.add_argument("-t","--threads", type=int, default=4)
    aps.add_argument("--stat", choices=["bp","bins"], default="bins",
                     help="Statistical unit for Fisher's exact test: base pairs (bp) or genomic bins (bins). Default: bins.")
    aps.add_argument("--score", choices=["giggle","odds"], default="giggle",
                     help="Ranking score: -log10(p)*log2(OR) (giggle) or -log10(p)*OR (odds).")
    aps.add_argument("--bin-size", type=int, default=4096,
                     help="Bin size in bp for --stat bins and for cache keys. Default: 4096.")
    aps.add_argument("--normalize-chr", action="store_true",
                     help="Normalize chromosome names in the QUERY to match the genome file (chr vs no-chr).")
    aps.add_argument("--cache-bins", action="store_true",
                     help="Load cached per-file bin masks if present, and save newly computed ones.")
    aps.add_argument("--fisher", choices=["one","two"], default="two",
                     help="Fisher's exact test: one-sided ('one'=greater) or two-sided ('two'). Default: two.")
    aps.add_argument("--fdr", action="store_true",
                     help="Calculate Benjamini–Hochberg FDR (Qvalue column).")

    args = ap.parse_args()
    if args.cmd == "index":
        index_db(args.db_dir, args.genome, threads=args.threads,
                 precompute_bins=args.precompute_bins, pre_bin_size=args.bin_size)
    elif args.cmd == "search":
        search(args.query_bed, args.db_dir, args.genome, out_tsv=args.out,
               threads=args.threads, stat=args.stat, score_mode=args.score,
               bin_size=args.bin_size, normalize_chr=args.normalize_chr,
               cache_bins=args.cache_bins, fisher_mode=args.fisher, do_fdr=args.fdr)

if __name__ == "__main__":
    main()
