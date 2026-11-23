Racial and ethnic representation in trials cited by the 2025 US hypertension and 2026 US dyslipidemia guidelines
================================================================================

Data and code to reproduce the numerical analyses in the research letter, including the values in both figure panels.

Run
---
Requires Python 3 and NumPy.

    python CODE-01-reproduce.py    prints reporting counts, participation ratios, and sensitivity analyses

Files
-----
CODE-01-reproduce.py          Reporting counts and pooled participation ratios, with all sensitivity analyses.
DATA-01-trials.csv            The 81 randomized trials newly cited by either guideline, one row per trial publication.
DATA-02-census-benchmarks.csv Census shares of US adults by group, in percent, at ages 20+ and 40+.
DATA-03-category-mismatches.csv  Trial categories that differ from the Census category they are compared with
                              (used only for the matching-categories sensitivity analysis).

DATA-01 columns
---------------
trial, doi, nct, other_registry   Trial name, publication DOI, and registry identifiers.
guidelines                    Citing guideline: HTN, LIPID, or both.
year                          Publication year.
geo_class                     US-only, US-inclusive multinational, non-US single-country, or non-US multinational.
us_enrolling                  True if at least one enrolling site was in the United States.
race_reported, hispanic_reported  True if the trial reported any explicit race or Hispanic/Latinx information.
race_denom, hisp_denom        Denominators for race counts and for the Hispanic/Latinx count, as reported.
A_WHITE ... A_HISPANIC        Participants reported in each group. Blank means not reported; 0 means reported as zero.
pooling_exclusion             Reason a trial's counts cannot be pooled (blank if usable).
note                          Population or category exceptions, such as ethnicity reported for a subset of participants.

Methods in brief
----------------
Counts come from ClinicalTrials.gov baseline tables (archived November 2026) and the trial publications.
Pooled analyses use trials with a US site and usable counts. For each group, the logit of the enrolled
proportion (0.5 added to each cell) is pooled with DerSimonian-Laird random effects and Hartung-Knapp 95%
confidence intervals, then divided by the group's share of US adults. Benchmarks are US Census Bureau
Vintage 2025 national estimates for July 1, 2025 (file nc-est2025-n-file12): race alone of any
ethnicity, and Hispanic/Latinx of any race. Counts from multinational trials include participants enrolled
outside the United States.
