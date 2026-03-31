#!/usr/bin/env python3
"""
CYP3A5 Tacrolimus Dosing — CPIC Guidelines
===========================================
Windows-compatible. Pure Python standard library only (NO cyvcf2 needed).
Handles both 'chr7' and '7' chromosome formats automatically.

Just run:
    python cyp3a5_tacrolimus.py

CPIC SNPs (GRCh37):
  *3  rs776746    chr7:99270539  G>A
  *6  rs10264272  chr7:99245974  C>T
  *7  rs41303343  chr7:99245899  T>TA (insertion)

References:
  Birdwell et al. Clin Pharmacol Ther 2015;98(1):19-24
  https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-cyp3a5/
"""

import csv
import glob
import gzip
import sys

# ── CPIC variant definitions ────────────────────────────────────────────────
VARIANTS = {
    "rs776746":   {"star": "*3", "chrom": "7", "pos": 99270539, "ref": "G", "alt": "A"},
    "rs10264272": {"star": "*6", "chrom": "7", "pos": 99245974, "ref": "C", "alt": "T"},
    "rs41303343": {"star": "*7", "chrom": "7", "pos": 99245899, "ref": "T", "alt": "TA"},
}

# Position lookup: normalised chrom + pos → rsid
POS_LOOKUP = {(v["chrom"], v["pos"]): rsid for rsid, v in VARIANTS.items()}


def norm_chrom(c):
    """Strip 'chr' prefix so 'chr7' == '7'."""
    return c.lower().replace("chr", "")


def find_vcf():
    for pattern in ("*.vcf.gz", "*.vcf", "*.bcf"):
        hits = glob.glob(pattern)
        if hits:
            print(f"[INFO] Found VCF: {hits[0]}")
            return hits[0]
    sys.exit("[ERROR] No VCF/VCF.gz file found in the current directory.")


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def resolve_base(allele_idx, ref, alts):
    if allele_idx == 0:
        return ref
    idx = allele_idx - 1
    return alts[idx] if idx < len(alts) else "?"


def classify(alleles):
    """
    alleles: list of 2 star strings e.g. ['*1','*3']
    Returns (phenotype, metabolizer_class, dose_recommendation)
    """
    star1 = alleles.count("*1")
    non_wt = [a for a in alleles if a != "*1"]

    if star1 == 2:
        pheno  = "Extensive Metabolizer"
        mclass = "Expresser — *1/*1"
    elif star1 == 1:
        pheno  = "Intermediate Metabolizer"
        mclass = f"Expresser — *1/{non_wt[0]}"
    else:
        pheno  = "Poor Metabolizer"
        mclass = f"Non-expresser — {alleles[0]}/{alleles[1]}"

    if star1 >= 1:
        dose = (
            "INCREASE starting dose to 1.5-2x the standard dose "
            "(max 0.3 mg/kg/day). Titrate to trough target with frequent TDM."
        )
    else:
        dose = (
            "Use STANDARD starting dose per transplant protocol "
            "(e.g. 0.1-0.15 mg/kg/day for kidney Tx). Titrate with routine TDM."
        )

    return pheno, mclass, dose


def main():
    vcf_path = find_vcf()
    samples = []
    sample_data = {}
    found = set()

    with open_vcf(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("#CHROM"):
                cols = line.split("\t")
                samples = cols[9:]
                sample_data = {s: {r: [] for r in VARIANTS} for s in samples}
                print(f"[INFO] {len(samples)} sample(s): {samples}")
                continue

            if line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            chrom_raw  = cols[0]
            pos        = int(cols[1])
            rsid_vcf   = cols[2]
            ref        = cols[3]
            alts       = cols[4].split(",")
            format_    = cols[8]
            chrom_norm = norm_chrom(chrom_raw)

            # Match by rsID first, then position
            matched = None
            if rsid_vcf in VARIANTS:
                matched = rsid_vcf
            else:
                matched = POS_LOOKUP.get((chrom_norm, pos))

            if matched is None:
                continue

            found.add(matched)
            print(f"[DEBUG] Matched {matched} ({VARIANTS[matched]['star']}) "
                  f"at {chrom_raw}:{pos}  REF={ref}  ALT={cols[4]}")

            fmt_fields = format_.split(":")
            gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else 0

            for i, sample in enumerate(samples):
                if 9 + i >= len(cols):
                    continue
                gt_raw = cols[9 + i].split(":")[gt_idx].replace("|", "/")
                bases = []
                for p in gt_raw.split("/"):
                    if p in (".", ""):
                        bases.append("?")
                    else:
                        bases.append(resolve_base(int(p), ref, alts))
                sample_data[sample][matched] = bases

    missing = set(VARIANTS) - found
    if missing:
        names = ", ".join(f"{r}({VARIANTS[r]['star']})" for r in missing)
        print(f"[WARN] Not found in VCF (assumed WT/*1): {names}")

    rows = []
    for sample in samples:
        no_func = []
        detail_parts = []

        for rsid, vdef in VARIANTS.items():
            bases = sample_data[sample][rsid]
            star  = vdef["star"]
            ref_b = vdef["ref"]

            if not bases:
                detail_parts.append(f"{star}({rsid}):not_found(WT assumed)")
                continue

            alt_count = sum(1 for b in bases if b != ref_b and b != "?")
            no_func.extend([star] * alt_count)
            detail_parts.append(f"{star}({rsid}):{'/'.join(bases)}")

        no_func   = no_func[:2]
        wt_count  = max(0, 2 - len(no_func))
        alleles   = sorted(["*1"] * wt_count + no_func)
        alleles   = (alleles + ["*1", "*1"])[:2]

        diplotype         = f"{alleles[0]}/{alleles[1]}"
        phenotype, mclass, dose = classify(alleles)

        rows.append({
            "sample":              sample,
            "diplotype":           diplotype,
            "phenotype":           phenotype,
            "metabolizer_class":   mclass,
            "dose_recommendation": dose,
            "variant_detail":      " | ".join(detail_parts),
            "vcf_source":          vcf_path,
        })

    out = "cyp3a5_tacrolimus_report.csv"
    fields = ["sample", "diplotype", "phenotype", "metabolizer_class",
              "dose_recommendation", "variant_detail", "vcf_source"]

    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"\n[INFO] Report saved -> {out}")
    print(f"\n{'SAMPLE':<20} {'DIPLOTYPE':<12} {'PHENOTYPE':<28} DOSE")
    print("-" * 100)
    for r in rows:
        print(f"{r['sample']:<20} {r['diplotype']:<12} {r['phenotype']:<28} {r['dose_recommendation'][:55]}...")


if __name__ == "__main__":
    main()
