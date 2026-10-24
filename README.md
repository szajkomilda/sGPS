# sGPS - similar Genomic Pattern Search

sGPS is a Python-based tool for indexing and searching genomic interval datasets. It compares the intervals of a query file against a database using Fisher's exact test and produces ranked enrichment scores.
Conceptually, sGPS is related to GIGGLE, however, the indexing strategy differs: sGPS stores datasets as bgzip+tabix files with an explicit manifest and optional precomputed bin masks (rather than GIGGLE's unified, B+ tree-based index). In addition, the statistical model uses Fisher's exact test (bin or basepair modes) with alternative scoring metrics and optional Benjamini-Hochberg FDR correction.
This repository accompanies a research manuscript and is intended for broader scientific use.

## Installation
Requirements: Python ≥ 3.8, pysam, numpy, scipy

Dependencies can be installed via pip or conda. The script runs directly and does not require compilation.

```
pip install pysam numpy scipy
```

## Usage overview
sGPS operates in two modes: index and search.
### 1. Index mode
Prepare a directory of BED/BED.GZ files for fast searching.

Usage:
python sGPS.py index <DB_DIR> --genome <FILE> [OPTIONS]

Options:

```
--genome <FILE> (required): Genome size file (chrom<TAB>length).
--threads <INT> (default: 4): Number of worker threads.
--precompute-bins (default: off): Precomputes per-file bin masks. Use the same --bin-size and --cache-bins when searching.
--bin-size <INT> (default: 4096): Bin size in bp for precomputation. 

```

Artifacts produced:
-bgzip+tabix-indexed BED files
-Manifest file (sGPS_manifest.json)
-Optionally cached bin masks (if --precompute-bins is used)

The .bed.gz and .tbi files are written next to the original BED files. Run index on a clean directory.

### 2. Search mode
Search a query BED file against an indexed database

Usage:
python sGPS.py search <query.bed> <DB_DIR> --genome <FILE> > <results.tsv>

Options:
```
--genome <FILE> (required): Genome size file (build must match query + DB).
--out <FILE> (default: "-" → stdout): Output TSV.
--threads <INT> (default: 4): Number of worker threads.
--stat {bins,bp} (default: bins): Statistical mode, genome bins or basepair coverage.
--score {giggle,odds} (default: giggle): 
  Ranking score:
    giggle: -log10(pval) * log2(odds_ratio), as in GIGGLE
    odds: -log10(pval) * odds_ratio
--bin-size <INT> (default: 4096): Bin size (bp) for bin mode.
--normalize-chr (default: off): Adds/strips the chr prefix in the QUERY to match the genome file (database files are not changed).
--cache-bins (default: off): Loads cached bin masks if present and saves newly computed ones.
--fisher {one,two} (default: two): Fisher test direction, one-sided (greater) or two-sided.
--fdr (default: off): Computes BH-FDR (q-values).
```

## Input data

### Genome file
Two-column tab-separated text file listing chromosome names and lengths.
Genome build is not fixed, but must match between query and database.

### BED files
Standard BED3+ format. BEDs may be bgzip-compressed. Query and database must use the same genome build. 

## Output
The search produces a ranked tab-separated table sorted by score (descending).

For bin-based mode:
```
  DB_File: Database file path
  OverlapBins: number of bins containing both query and database intervals
  QueryIntervals: number of intervals in the query
  DBIntervals: number of intervals in the database file
  OddsRatio: Fisher's exact test odds ratio
  Pvalue: Fisher's exact test p-value
  ComboScore: ranking score; with --score giggle it uses the GIGGLE formula, but computed from bin- or bp-based tables, so values are not comparable to GIGGLE output.
  Qvalue: Benjamini-Hochberg FDR (if --fdr flag used)
```

For basepair mode:
Analogous columns with bp-based coverage metrics instead of bin counts.

## Reproducibility Notes
Minimal requirements for reproducibility:
```
  Query BED file
  BED database
  Genome size file
```
```
   python sGPS.py index data/db --genome data/genome.txt --precompute-bins
   python sGPS.py search data/query.bed data/db --genome data/genome.txt --cache-bins --fdr -o results.tsv
```

## License
Distributed under the MIT License. See LICENSE file for details.

## Citation
If you use sGPS in scientific work, please cite the associated manuscript.

GIGGLE: Layer, R. M., Pedersen, B. S., DiSera, T., Marth, G. T., Gertz, J., & Quinlan, A. R. (2018). GIGGLE: a search engine for large-scale integrated genome analysis. Nature methods, 15(2), 123–126. https://doi.org/10.1038/nmeth.4556

## Contact
For questions, issues, or contributions, please open an issue on GitHub.

