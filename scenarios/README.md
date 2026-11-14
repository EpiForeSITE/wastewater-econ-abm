# Scenario comparison


This page compares the scenarios in [scenarios.yml](../scenarios.yml),
using the outputs in [results/](../results). Run `make` from the
repository root to render every scenario and this page. Each scenario
has its own report:

| scenario | response | description |
|:---|:---|:---|
| [baseline](baseline.md) | cases | Case-based surveillance: reported cases start a public health alert and school exclusions; wastewater samples are recorded but not acted on. |
| [wastewater](wastewater.md) | wastewater | Baseline plus wastewater surveillance as practiced in Utah: samples twice a week, with results four days after collection; positive results also start the alert. |
| [wastewater-fast](wastewater-fast.md) | wastewater | Wastewater surveillance with a faster turnaround: results two days after collection. |

All scenarios share 1000 replicates, the same population, vaccination
and calibrated transmission rate (0.115 per contact). Counts are for the
120,000-agent model; multiply by 10.17 for Salt Lake County.

## Outcomes by scenario

Means per replicate, with medians of infections; most introductions stay
small, so means are driven by the outbreaks that grow.

| scenario | infections_mean | infections_median | infections_q90 | admissions | hospital_days | excluded_student_days | samples | alert_days | false_positives |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 32.6 | 4 | 95.1 | 0.6 | 3.9 | 13930.1 | 105 | 88.5 | 0.5 |
| wastewater | 32.2 | 4 | 95.1 | 0.6 | 4.0 | 13970.1 | 105 | 107.5 | 0.5 |
| wastewater-fast | 31.8 | 4 | 94.1 | 0.6 | 3.9 | 13928.0 | 105 | 107.6 | 0.5 |

## Paired differences from the baseline

Replicate *i* of every scenario starts from the same school and the same
random draws, so differences within a replicate are effects of the
scenario’s settings. “Averted” is baseline minus scenario; positive
values favor the scenario.

| scenario | infections_averted | infections_averted_median | fewer_infections_pct | more_infections_pct | admissions_averted | hospital_days_averted | excluded_days_averted | alert_days_added |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| wastewater | 0.4 | 0 | 4.9 | 3.3 | 0 | 0.0 | -40.0 | 19.0 |
| wastewater-fast | 0.8 | 0 | 5.7 | 3.8 | 0 | 0.1 | 2.1 | 19.2 |

## Timeliness

For replicates with at least one reported case, the lead time is the day
of the first case report minus the day the first true-positive
wastewater result arrived (positive values mean wastewater was first).
False positives are excluded.

| scenario | replicates | ww_positive_pct | ww_first_pct | lead_days_median | lead_days_q25 | lead_days_q75 |
|:---|---:|---:|---:|---:|---:|---:|
| baseline | 704 | 62.6 | 21.9 | -7 | -31 | 5 |
| wastewater | 716 | 62.3 | 22.6 | -6 | -30 | 5 |
| wastewater-fast | 721 | 62.7 | 28.0 | -4 | -28 | 5 |

## Alert timing

Every scenario raises case detection during a public health alert; the
wastewater scenarios can start it before the first case report. Days
gained compares, within each replicate, the day the alert first started
in the baseline (first case report) and in the scenario. Replicates in
which the scenario’s first alert followed a false positive before the
outbreak was detected also count as alerts started by wastewater.

| scenario | replicates | started_by_ww_pct | earlier_pct | days_gained_q50 | days_gained_q75 | days_gained_q90 |
|:---|---:|---:|---:|---:|---:|---:|
| wastewater | 704 | 25.9 | 25.9 | 0 | 1 | 16 |
| wastewater-fast | 704 | 30.7 | 30.7 | 0 | 3 | 18 |
