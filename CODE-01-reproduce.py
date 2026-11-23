"""Reproduce every number reported in the research letter.

Run:  python CODE-01-reproduce.py
Reads DATA-01 to DATA-03 from this folder and prints the results in the order they appear in the letter.
"""
import csv
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GROUPS = ["WHITE", "BLACK", "ASIAN", "AIAN", "HISPANIC"]


def read(name):
    with open(HERE / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


TRIALS = read("DATA-01-trials.csv")
BENCH = {r["group"]: r for r in read("DATA-02-census-benchmarks.csv")}
# trial and group pairs whose published category differs from the Census category it is compared with
MISMATCH = {(r["nct"], r["group"]) for r in read("DATA-03-category-mismatches.csv")
            if r["census_benchmark"] != r["group"]}
US_SITE = [r for r in TRIALS if r["us_enrolling"] == "True"]
POOLED = [r for r in US_SITE if not r["pooling_exclusion"] and any(r["A_" + g] for g in GROUPS)]


def t_critical(df):
    """Two-sided 95% Student t quantile, by Gauss-Legendre integration of the density (numpy only)."""
    nodes, weights = np.polynomial.legendre.leggauss(32)
    c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    lo, hi = 0.0, 8.0
    for _ in range(50):
        mid = (lo + hi) / 2
        x = (nodes + 1) * mid / 2
        area = c * mid / 2 * np.dot(weights, (1 + x * x / df) ** (-(df + 1) / 2))
        lo, hi = (mid, hi) if area < 0.475 else (lo, mid)
    return (lo + hi) / 2


def pool(trials, group, age="20"):
    """Random-effects pooled share (logit scale, 0.5 added to each cell, DerSimonian-Laird,
    Hartung-Knapp 95% CI) divided by the Census share of US adults."""
    denom = "hisp_denom" if group == "HISPANIC" else "race_denom"
    x, n = np.array([(float(r["A_" + group]), float(r[denom])) for r in trials
                     if r["A_" + group] and r[denom]]).T
    y = np.log((x + 0.5) / (n - x + 0.5))
    v = 1 / (x + 0.5) + 1 / (n - x + 0.5)
    w = 1 / v
    q = np.sum(w * (y - np.average(y, weights=w)) ** 2)
    tau2 = max(0.0, (q - (len(x) - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    w = 1 / (v + tau2)
    mu = np.average(y, weights=w)
    half = t_critical(len(x) - 1) * math.sqrt(np.average((y - mu) ** 2, weights=w) / (len(x) - 1))
    pct = 100 / (1 + np.exp(-np.array([mu, mu - half, mu + half])))
    bench = float(BENCH[group]["pct_us_adults_%s_plus" % age])
    return dict(group=group, trials=len(x), zero_cells=int((x == 0).sum()), events=int(x.sum()),
                participants=int(n.sum()), pooled_pct=pct[0], benchmark_pct=bench,
                ratio=pct[0] / bench, lower=pct[1] / bench, upper=pct[2] / bench,
                combined_ratio=100 * x.sum() / n.sum() / bench)


def ci(f):
    return "%.2f (95%% CI, %.2f-%.2f)" % (f["ratio"], f["lower"], f["upper"])


def report():
    def reporting(label, rows):
        race = sum(r["race_reported"] == "True" for r in rows)
        eth = sum(r["hispanic_reported"] == "True" for r in rows)
        print("%-22s n=%-3d race %d (%.0f%%), ethnicity %d (%.0f%%)"
              % (label, len(rows), race, 100 * race / len(rows), eth, 100 * eth / len(rows)))

    print("Reporting (Figure A)")
    reporting("All trials", TRIALS)
    reporting("With a US site", US_SITE)
    reporting("Without a US site", [r for r in TRIALS if r["us_enrolling"] == "False"])

    print("\nPooled trials: %d of %d with a US site" % (len(POOLED), len(US_SITE)))
    for r in US_SITE:
        if r not in POOLED:
            print("  not pooled:", r["trial"], "-", r["pooling_exclusion"] or "no race or ethnicity counts")

    print("\nParticipation ratios, Census adults 20 years or older (Figure B)")
    for g in GROUPS:
        f = pool(POOLED, g)
        print("  %-8s trials %2d  pooled %4.1f%% vs %4.1f%%  ratio %s" % (g, f["trials"], f["pooled_pct"], f["benchmark_pct"], ci(f)))
    h = pool(POOLED, "HISPANIC")
    print("  Hispanic/Latinx participants combined across trials: %s of %s (%.1f%%), ratio %.2f"
          % (format(h["events"], ","), format(h["participants"], ","), 100 * h["events"] / h["participants"], h["combined_ratio"]))
    a = pool(POOLED, "AIAN")
    print("  American Indian or Alaska Native: %d of %d trials reported zero participants" % (a["zero_cells"], a["trials"]))

    print("\nSensitivity analyses")
    print("  Hispanic/Latinx, Census adults 40 years or older:", ci(pool(POOLED, "HISPANIC", "40")))
    for g in ["HISPANIC", "ASIAN", "AIAN"]:
        print("  %-8s matching categories only: %s" % (g, ci(pool([r for r in POOLED if (r["nct"], g) not in MISMATCH], g))))

    print("\nBlack participants by enrollment geography")
    for label, rows in [("US-only", [r for r in POOLED if r["geo_class"] == "US-only"]),
                        ("Multinational with a US site", [r for r in POOLED if r["geo_class"] != "US-only"])]:
        f = pool(rows, "BLACK")
        print("  %-29s trials %2d  ratio %s" % (label, f["trials"], ci(f)))


if __name__ == "__main__":
    report()
