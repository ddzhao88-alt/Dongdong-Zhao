#!/usr/bin/env python3
from __future__ import annotations

import argparse, csv, gzip, os, re
from typing import Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import xlsxwriter

# -------- Parameters (match the manuscript analysis logic) --------
POS_RANGE  = range(43, 58)              # start positions 43–57
DIST_RANGE = [p - 1 for p in POS_RANGE] # offsets 42–56 after spacer end
TN_LEN = 20

GS_UPSTREAM_START = 25   # -25
GS_UPSTREAM_END   = 5    # -5
GS_LEN = GS_UPSTREAM_START - GS_UPSTREAM_END  # 20


# -------- Helpers --------
def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement()).upper()

def open_fq(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def norm(col: str) -> str:
    return re.sub(r"\W+", "_", col.lstrip("\ufeff").strip().lower())

def ed_leq1_distance(a: str, b: str) -> int:
    """Return edit distance if <=1 (sub/ins/del), else return 2."""
    a = a.upper(); b = b.upper()
    la, lb = len(a), len(b)
    if abs(la - lb) > 1:
        return 2

    # same length: <=1 substitution
    if la == lb:
        mm = 0
        for x, y in zip(a, b):
            if x != y:
                mm += 1
                if mm > 1:
                    return 2
        return mm

    # length differs by 1: single insertion/deletion
    longer, shorter = (a, b) if la > lb else (b, a)
    i = j = 0
    edits = 0
    while i < len(longer) and j < len(shorter):
        if longer[i] == shorter[j]:
            i += 1; j += 1
        else:
            edits += 1
            if edits > 1:
                return 2
            i += 1
    return 1  # <=1 indel

def ed_leq1(a: str, b: str) -> bool:
    return ed_leq1_distance(a, b) <= 1

def seed_positions(L: int, seed_len: int) -> list[int]:
    if L <= seed_len:
        return [0]
    pos = {0, (L - seed_len)//2, L - seed_len}
    pos.add(max(0, (L - seed_len)//3))
    pos.add(max(0, 2*(L - seed_len)//3))
    return sorted(pos)

def approx_find_ed1(seq: str, pattern: str) -> Tuple[int, int]:
    """Find pattern in seq with ED<=1. Return (best_idx, best_ed) or (-1,2)."""
    s = seq.upper()
    pat = pattern.upper()
    L = len(pat)
    if L == 0 or len(s) < L - 1:
        return -1, 2

    idx_exact = s.find(pat)
    if idx_exact != -1:
        return idx_exact, 0

    seed_len = 10 if L >= 30 else 8
    seed_len = min(seed_len, L)

    cands = set()
    for p0 in seed_positions(L, seed_len):
        seed = pat[p0:p0+seed_len]
        start = 0
        while True:
            hit = s.find(seed, start)
            if hit == -1:
                break
            est = hit - p0
            for st in (est - 1, est, est + 1):
                if st < 0:
                    continue
                if st + (L - 1) <= len(s):
                    cands.add(st)
            start = hit + 1

    best_idx = -1
    best_ed = 2

    for st in sorted(cands):
        for wlen in (L, L - 1, L + 1):
            if wlen <= 0 or st + wlen > len(s):
                continue
            window = s[st:st+wlen]
            ed = ed_leq1_distance(window, pat)
            if ed <= 1:
                if ed < best_ed or (ed == best_ed and (best_idx == -1 or st < best_idx)):
                    best_ed = ed
                    best_idx = st
                    if best_ed == 0:
                        return best_idx, best_ed

    # brute fallback
    max_start = len(s) - (L - 1)
    for st in range(0, max_start + 1):
        for wlen in (L, L - 1, L + 1):
            if wlen <= 0 or st + wlen > len(s):
                continue
            window = s[st:st+wlen]
            ed = ed_leq1_distance(window, pat)
            if ed <= 1:
                if ed < best_ed or (ed == best_ed and (best_idx == -1 or st < best_idx)):
                    best_ed = ed
                    best_idx = st
                    if best_ed == 0:
                        return best_idx, best_ed

    return best_idx, best_ed

def approx_contains_ed1(seq: str, pattern: str) -> bool:
    idx, ed = approx_find_ed1(seq, pattern)
    return idx != -1 and ed <= 1

def orient_by_spacer(seq: str, spacer: str) -> Tuple[str, int]:
    """Orient read by spacer with ED<=1, return (oriented_seq, spacer_idx) or ("",-1)."""
    s = seq.upper()
    sp = spacer.upper()

    idx, ed = approx_find_ed1(s, sp)
    if idx != -1 and ed <= 1:
        return s, idx

    s_rc = revcomp(s)
    idx2, ed2 = approx_find_ed1(s_rc, sp)
    if idx2 != -1 and ed2 <= 1:
        return s_rc, idx2

    return "", -1

def derive_genomic_site(genome_seq: str, spacer: str) -> str:
    """From genome_seq, locate spacer (ED<=1; allow RC) and return -25~-5 (20bp)."""
    g = genome_seq.upper()
    sp = spacer.upper()

    idx, ed = approx_find_ed1(g, sp)
    if idx == -1:
        g = revcomp(g)
        idx, ed = approx_find_ed1(g, sp)
        if idx == -1:
            raise ValueError("Cannot find spacer in genome_seq (both strands, ED<=1).")

    if idx < GS_UPSTREAM_START:
        raise ValueError(f"Spacer upstream <{GS_UPSTREAM_START} bp; cannot derive -25~-5 window.")
    gs = g[idx - GS_UPSTREAM_START : idx - GS_UPSTREAM_END]
    if len(gs) != GS_LEN:
        raise ValueError("Derived genomic_site length != 20 bp.")
    return gs

def tally(fq_path: str, tn20: str, spacer: str, genomic_site: str):
    tn = tn20.upper()
    sp = spacer.upper()
    gs = genomic_site.upper()

    total = 0
    oriented = 0
    integ = 0
    unedited = 0
    waste = 0
    dist = {d: 0 for d in DIST_RANGE}

    with open_fq(fq_path) as fh:
        for rec in SeqIO.parse(fh, "fastq"):
            total += 1
            raw_seq = str(rec.seq)

            seq, idx = orient_by_spacer(raw_seq, sp)
            if idx == -1:
                waste += 1
                continue
            oriented += 1

            # integrated if donor tag found at offsets 42–56 after spacer end (allow <=1 edit; allow +/-1 indel length)
            found = False
            for d in DIST_RANGE:
                p = idx + len(sp) + d
                for wlen in (TN_LEN, TN_LEN - 1, TN_LEN + 1):
                    if wlen <= 0 or p + wlen > len(seq):
                        continue
                    window = seq[p:p+wlen]
                    if ed_leq1(window, tn):
                        found = True
                        integ += 1
                        dist[d] += 1
                        break
                if found:
                    break
            if found:
                continue

            # unedited if genomic_site present (ED<=1)
            if approx_contains_ed1(seq, gs):
                unedited += 1

    return total, oriented, integ, unedited, waste, dist

def read_single_sample(csv_path: str) -> Dict[str, str]:
    with open(csv_path, newline="") as f:
        rdr = csv.DictReader(f)
        if rdr.fieldnames is None:
            raise ValueError("sample_info.csv has no header.")
        nf = {c: norm(c) for c in rdr.fieldnames}

        def get_key(*cands):
            for k, v in nf.items():
                if v in cands:
                    return k
            raise ValueError(f"Missing required column: {cands}")

        name_k = get_key("name", "sample", "sample_name")
        re_k   = get_key("donor_re", "re", "donorrightend")
        cr_k   = get_key("crrna", "spacer", "crrna_spacer", "guide", "grna")
        gn_k   = get_key("genome_seq", "genomic_seq", "genomic", "genome", "target_seq")

        rows = list(rdr)
        if len(rows) == 0:
            raise ValueError("sample_info.csv contains no rows.")
        if len(rows) > 1:
            raise ValueError("This demo script expects exactly 1 sample in sample_info.csv.")

        row = rows[0]
        return {
            "name": row[name_k].strip(),
            "donor_re": row[re_k].strip().upper(),
            "spacer": row[cr_k].strip().upper(),
            "genome_seq": row[gn_k].strip().upper(),
        }

def main():
    ap = argparse.ArgumentParser(description="MaCAST integration tally (demo-friendly).")
    ap.add_argument("--csv", required=True, help="sample_info.csv (single sample)")
    ap.add_argument("--fastq", required=True, help="FASTQ(.gz) for this sample")
    ap.add_argument("--excel", default="output.xlsx", help="output .xlsx (default: output.xlsx)")
    args = ap.parse_args()

    if not os.path.isfile(args.fastq):
        raise FileNotFoundError(f"FASTQ not found: {args.fastq}")

    info = read_single_sample(args.csv)
    if len(info["donor_re"]) < TN_LEN:
        raise ValueError("donor_re length < 20 bp; cannot derive donor tag.")

    tn20 = info["donor_re"][:TN_LEN]
    genomic_site = derive_genomic_site(info["genome_seq"], info["spacer"])

    total, oriented, integ, unedited, waste, dist = tally(
        fq_path=args.fastq,
        tn20=tn20,
        spacer=info["spacer"],
        genomic_site=genomic_site,
    )

    denom = integ + unedited
    eff = 100.0 * integ / denom if denom else 0.0

    wb = xlsxwriter.Workbook(args.excel)
    ws = wb.add_worksheet("Summary")
    bold = wb.add_format({"bold": True})

    header = [
        "Name","Total_reads","Oriented_reads","Integrated_reads",
        "Unedited_reads","Waste_reads","Integration_efficiency (%)"
    ] + [f"{d} bp" for d in DIST_RANGE]
    ws.write_row(0, 0, header, bold)

    ws.write_row(
        1, 0,
        [info["name"], total, oriented, integ, unedited, waste, f"{eff:.2f}"]
        + [dist[d] for d in DIST_RANGE]
    )

    wb.close()
    print(f"Done. Wrote: {args.excel}")

if __name__ == "__main__":
    main()
