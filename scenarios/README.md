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

| scenario | infections_mean | infections_median | infections_q90 | admissions | hospital_days | excluded_student_days | samples | alert_days | false_positives | campaign_doses |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 9.5 | 4 | 25.0 | 0.2 | 1.2 | 3675.1 | 105 | 54.1 | 0.6 | 5405.8 |
| wastewater | 8.4 | 4 | 22.1 | 0.2 | 1.2 | 3396.3 | 105 | 74.2 | 0.6 | 6848.8 |
| wastewater-fast | 8.3 | 4 | 21.0 | 0.2 | 1.1 | 3373.7 | 105 | 75.3 | 0.6 | 6881.9 |

## Paired differences from the baseline

Replicate *i* of every scenario starts from the same school and the same
random draws, so differences within a replicate are effects of the
scenario’s settings. “Averted” is baseline minus scenario; positive
values favor the scenario.

| scenario | infections_averted | infections_averted_median | fewer_infections_pct | more_infections_pct | admissions_averted | hospital_days_averted | excluded_days_averted | alert_days_added | doses_added |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| wastewater | 1.1 | 0 | 10.6 | 3.5 | 0 | 0.1 | 278.8 | 20.1 | 1443.1 |
| wastewater-fast | 1.3 | 0 | 12.3 | 3.9 | 0 | 0.2 | 301.4 | 21.2 | 1476.1 |

## Timeliness

For replicates with at least one reported case, the lead time is the day
of the first case report minus the day the first true-positive
wastewater result arrived (positive values mean wastewater was first).
False positives are excluded.

| scenario | replicates | ww_positive_pct | ww_first_pct | lead_days_median | lead_days_q25 | lead_days_q75 |
|:---|---:|---:|---:|---:|---:|---:|
| baseline | 704 | 56.5 | 21.9 | -3.0 | -18.8 | 8 |
| wastewater | 714 | 56.0 | 22.4 | -2.5 | -18.0 | 5 |
| wastewater-fast | 722 | 56.5 | 28.1 | 0.0 | -16.0 | 7 |

## Alert timing

Every scenario raises case detection and runs a vaccination campaign
during a public health alert; the wastewater scenarios can start it
before the first case report. Days gained compares, within each
replicate, the day the alert first started in the baseline (first case
report) and in the scenario. Replicates in which the scenario’s first
alert followed a false positive before the outbreak was detected also
count as alerts started by wastewater.

| scenario | replicates | started_by_ww_pct | earlier_pct | days_gained_q50 | days_gained_q75 | days_gained_q90 |
|:---|---:|---:|---:|---:|---:|---:|
| wastewater | 704 | 25.9 | 25.9 | 0 | 1 | 16 |
| wastewater-fast | 704 | 30.7 | 30.7 | 0 | 3 | 18 |
