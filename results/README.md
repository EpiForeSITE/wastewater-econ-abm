# Simulation results

Each folder holds the outputs of one scenario in
[scenarios.yml](../scenarios.yml), written by `make <scenario>` (see the
[main README](../README.md)); do not edit them by hand. The rendered report
for each is in [scenarios/](../scenarios).

```
results/<scenario>/
  replicates.csv         one row per replicate: outcomes for the evaluation
  daily.csv.gz           one row per replicate and day
  infections_by_age.csv  one row per replicate and age band
  metadata.csv           settings, calibration, and commits
```

## Study design

Each replicate seeds one infection in a student of a school drawn uniformly
at random, then runs 365 days. Replicate *i* uses the same seed, origin
school, population and vaccination in every scenario, so scenarios are
identical until their settings make them diverge. Compare scenarios
**within** replicates (paired differences, joined on `replicate`), as in
[scenarios/README.md](../scenarios/README.md).

## Units and scaling

All counts are for the 120,000-agent model. Multiply counts by
`county_scale` in `metadata.csv` (about 10.17) to approximate Salt Lake
County totals. Wastewater samples, alerts and alert days do not scale: there
is one catchment for the county.

Most introductions die out, so outcomes are very skewed. Report medians and
upper quantiles alongside means.

## `replicates.csv`

| Column | Meaning |
|:---|:---|
| `scenario` | Scenario name in `scenarios.yml` |
| `response` | `cases` (alerts and exclusions follow reported cases) or `wastewater` (positive samples can also start alerts) |
| `replicate` | Replicate number; the key for paired comparisons |
| `seed` | Random seed of the replicate |
| `seed_school` | Index of the origin school (0-based) |
| `infections` | Agents ever infected, including the seed |
| `rash_cases` | Rash onsets by day 365 |
| `detected_cases` | Cases detected (and isolated) |
| `reported_cases` | Detected cases reported to public health |
| `first_report_day` | Day of the first case report (empty if none) |
| `admissions` | Hospital admissions |
| `hospital_days` | Hospital bed-days |
| `isolation_days` | Person-days detected cases spend isolated at home |
| `last_infection_day` | Day of the last new infection (0 if none) |
| `active_at_end` | Whether any agent is still infected on day 365 |
| `samples` | Wastewater samples taken |
| `positive_samples` | Positive samples, including false positives |
| `false_positive_samples` | Positive samples drawn as false positives |
| `alerts` | Distinct public health alerts (an alert extended by a new trigger counts once) |
| `first_alert_day` | Day the first alert started (empty if none) |
| `first_alert_trigger` | What started the first alert: `case` or `wastewater` (empty if none) |
| `first_true_positive_day` | Day the first true-positive wastewater result arrived: sample day plus lab delay (empty if none); recorded in every scenario |
| `alert_days` | Days with an alert on (case detection is higher) |
| `exclusion_episodes` | School exclusions started |
| `first_exclusion_day` | Day of the first exclusion (empty if none) |
| `excluded_student_days` | Student-days excluded from school |

An alert raises the daily detection probability so that
`alert_reporting_fraction` (50%) of cases are reported instead of 25%. With
the `cases` response only reported cases start alerts; wastewater samples
are still recorded.

## `daily.csv.gz`

| Column | Meaning |
|:---|:---|
| `scenario`, `replicate` | Keys, as above |
| `day` | Day 1 to 365 |
| `new_infections` | New infections that day |
| `new_rash` | New rash onsets |
| `new_detected` | Cases detected and isolated |
| `reported` | Cases reported (they trigger school exclusions) |
| `admissions` | New hospital admissions |
| `hospitalized` | Agents in hospital at the end of the day |
| `isolated` | Detected cases isolated at home at the end of the day |
| `excluded` | Students excluded from school |
| `alert_on` | Whether a public health alert is on |
| `ww_obs` | Measured concentration (copies/mL); empty on non-sampling days |
| `ww_positive` | Whether the sample was positive; empty on non-sampling days |
| `ww_false` | Whether the positive was false; empty on non-sampling days |

`day` is the sampling day; the result becomes known `ww_lab_delay` days
later (see the scenario's report).

## `infections_by_age.csv`

`scenario`, `replicate`, `age_group` (five-year band, `0-4` to `75+`), and
`infections`, the agents ever infected in that band, including the seed.

## `metadata.csv`

`scenario`, `response`, `description`, `n` (agents), `county_population`,
`county_scale`, `replicates`, `schools`, `transmission_rate` (calibrated, per
contact), `simulated_r0`, `epiworldR_commit`, `git_commit`, and `rendered`
(time of the render).
