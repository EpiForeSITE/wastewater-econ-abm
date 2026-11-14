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
calibrated transmission rate (0.152 per contact). Counts are for the
120,000-agent model; multiply by 10.17 for Salt Lake County.

## Outcomes by scenario

Means per replicate, with medians of infections; most introductions stay
small, so means are driven by the outbreaks that grow.

| scenario | infections_mean | infections_median | infections_q90 | admissions | hospital_days | excluded_student_days | samples | alert_days | reduced_contact_days |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 38.2 | 4 | 140.1 | 6.1 | 43.2 | 24105.0 | 105 | 71.3 | 0 |
| wastewater | 16.2 | 4 | 45.1 | 2.5 | 17.2 | 13154.3 | 105 | 57.2 | 685860 |
| wastewater-fast | 16.3 | 4 | 46.6 | 2.7 | 18.2 | 13169.1 | 105 | 58.3 | 699660 |

## Paired differences from the baseline

Replicate *i* of every scenario starts from the same school and the same
random draws, so differences within a replicate are effects of the
scenario’s settings. “Averted” is baseline minus scenario; positive
values favor the scenario.

| scenario | infections_averted | infections_averted_median | fewer_infections_pct | more_infections_pct | admissions_averted | hospital_days_averted | excluded_days_averted | reduced_contact_days_added |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| wastewater | 22.0 | 0 | 24.5 | 14.5 | 3.6 | 26.0 | 10950.7 | 685860 |
| wastewater-fast | 21.9 | 0 | 25.5 | 15.5 | 3.5 | 24.9 | 10935.9 | 699660 |

## Timeliness

For replicates with at least one reported case, the lead time is the day
of the first case report minus the day of the first positive wastewater
result (positive values mean wastewater was first).

| scenario        | replicates | alert_first_pct | lead_days_median |
|:----------------|-----------:|----------------:|-----------------:|
| baseline        |        197 |             5.1 |              -27 |
| wastewater      |        197 |             5.1 |              -27 |
| wastewater-fast |        197 |             9.6 |              -25 |
