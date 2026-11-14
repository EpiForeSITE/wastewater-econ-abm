# Scenario comparison


This page compares the scenarios in [scenarios.yml](../scenarios.yml),
using the outputs in [results/](../results). Run `make` from the
repository root to render every scenario and this page. Each scenario
has its own report:

| scenario | response | description |
|:---|:---|:---|
| [baseline](baseline.md) | cases | School exclusions after reported student cases; no wastewater surveillance acts. |
| [wastewater](wastewater.md) | wastewater | Baseline plus wastewater surveillance as practiced in Utah: samples twice a week, with results four days after collection. |
| [wastewater-fast](wastewater-fast.md) | wastewater | Wastewater surveillance with a faster turnaround: results two days after collection. |

All scenarios share 200 replicates, the same population, vaccination and
calibrated transmission rate (0.115 per contact). Counts are for the
120,000-agent model; multiply by 10.17 for Salt Lake County.

## Outcomes by scenario

Means per replicate, with medians of infections; most introductions stay
small, so means are driven by the outbreaks that grow.

| scenario | infections_mean | infections_median | infections_q90 | admissions | hospital_days | excluded_student_days | samples | alert_days | reduced_contact_days |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 67.6 | 5 | 282.5 | 1.1 | 7.6 | 15072.4 | 105 | 92.7 | 0 |
| wastewater | 35.6 | 5 | 130.2 | 0.7 | 4.8 | 11530.9 | 105 | 83.3 | 999120 |
| wastewater-fast | 42.6 | 5 | 127.1 | 0.7 | 5.1 | 12544.6 | 105 | 87.6 | 1051740 |

## Paired differences from the baseline

Replicate *i* of every scenario starts from the same school and the same
random draws, so differences within a replicate are effects of the
scenario’s settings. “Averted” is baseline minus scenario; positive
values favor the scenario.

| scenario | infections_averted | infections_averted_median | fewer_infections_pct | more_infections_pct | admissions_averted | hospital_days_averted | excluded_days_averted | reduced_contact_days_added |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| wastewater | 32.1 | 0 | 27.5 | 16.0 | 0.4 | 2.9 | 3541.4 | 999120 |
| wastewater-fast | 25.0 | 0 | 25.0 | 20.5 | 0.4 | 2.6 | 2527.8 | 1051740 |

## Timeliness

For replicates with at least one reported case, the lead time is the day
of the first case report minus the day the first true-positive
wastewater result arrived (positive values mean wastewater was first).
False positives are excluded.

| scenario | replicates | ww_positive_pct | ww_first_pct | lead_days_median | lead_days_q25 | lead_days_q75 |
|:---|---:|---:|---:|---:|---:|---:|
| baseline | 141 | 67.4 | 19.9 | -9 | -32.5 | 3.5 |
| wastewater | 140 | 69.3 | 20.0 | -10 | -36.0 | 3.0 |
| wastewater-fast | 139 | 68.3 | 27.3 | -8 | -32.5 | 6.5 |
