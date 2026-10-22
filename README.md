# Measles wastewater surveillance with epiworldR


This project simulates measles outbreaks in a synthetic city of 10,000
people and compares three scenarios:

1.  **Baseline**: the outbreak runs with no surveillance-triggered
    intervention.
2.  **Wastewater (WW)**: wastewater is sampled three times a week,
    results are available after a lab turnaround time, and an
    intervention starts once a result reaches a threshold.
3.  **Cases**: each case is detected during its rash period with
    probability 70% and reported after a delay, and the same
    intervention starts once cumulative reported cases reach a
    threshold. Clinical detection of measles is far from perfect, so we
    vary this probability in a sensitivity analysis (see [Case detection
    probability](#case-detection-probability)).

The intervention reduces the non-household contacts (school, work, and
other settings) of school-age groups (5–19) by 80% for the rest of the
outbreak. The WW threshold is calibrated so both triggers respond to the
same number of cases; the difference between them is *when* each system
sees those cases.

## Model

### Population and contacts

- **Disease model**: `measles::ModelMeaslesMixing()` (latent, prodromal,
  rash, hospitalized, recovered), with one entity per 5-year age group
  (0–4, …, 75+). Individual-level isolation, quarantine, and contact
  tracing are off in every scenario, so the only response is the
  triggered intervention.
- **Age distribution**: Salt Lake County, UT (2025 Census estimates, via
  the
  [multigroup.vaccine](https://github.com/EpiForeSITE/multigroup-vaccine)
  package), rescaled to 10,000 people.
- **Contacts**: the setting-specific (household, school, work, other)
  contact matrices from
  [Epistorm-Mix](https://github.com/epistorm/Epistorm-Mix) (see
  [data/epistorm](data/epistorm/README.md)), re-weighted to the county’s
  age distribution and made reciprocal with
  `measles::make_cmat_symmetric()`.
- **Transmission**: the contact matrix is scaled to R0 = 7 (a
  parameter), first analytically and then by simulation, because the
  agent-based model realizes a lower R0 than the analytic formula. The
  simulated R0 is 6.91. Under the intervention it drops to about 5.4;
  school-age contacts outside the home account for only part of
  transmission.
- **Natural history**: 10 days latent, 4 days prodromal, and 4 days of
  rash. That puts rash onset about 14 days after infection and makes
  people infectious from 4 days before to 4 days after rash onset.
  People with rash keep 50% of their contacts. 20% are hospitalized.

### Wastewater shedding

The shedding model follows the structure of the CDC
[wwinference](https://github.com/CDCgov/ww-inference-model) model. The
expected concentration on day $t$ is

$$C(t) = \frac{1}{\alpha N} \sum_{\tau=0}^{\tau_{end}} s(\tau)\, G(t - \tau),$$

where $G(d)$ is the total number of genome copies shed by the people
whose prodrome started on day $d$, $\alpha$ is the wastewater volume per
person per day, and $N$ is the population. The shedding profile
$s(\tau)$ is piecewise linear on the log10 scale: it rises to its peak
at $\tau_{peak}$ and declines back to zero at $\tau_{end}$, normalized
to sum to one. Two changes from the COVID model:

- $\tau$ is counted from **prodrome onset**, not infection, because
  measles RNA appears in urine during the prodrome and peaks around rash
  onset.
- Each person’s total shedding is drawn from a log-normal distribution,
  since measured measles shedding varies over four orders of magnitude.

Measurements are log-normal around $C(t)$ and recorded as non-detects
below the limit of detection (LOD).

See [Parameters](#parameters) for the values and their sources. With
these values, one person with median shedding, on their peak day, in a
city of 10,000 produces a concentration near the LOD, roughly in line
with the STOTEN 2024 estimate that about 8 cases per 10,000 people are
needed for 50% detection.

``` r
k <- shedding_kernel(p$ww_t_peak, p$ww_t_end, p$ww_log10_vpeak)
ggplot(data.frame(tau = seq_along(k) - 1, s = k), aes(tau, s)) +
  geom_col(fill = "steelblue") +
  geom_vline(xintercept = p$prodromal_period, linetype = 2) +
  labs(x = "Days since prodrome onset", y = "Share of total shedding") +
  theme_minimal()
```

<div id="fig-kernel">

![](README_files/figure-commonmark/fig-kernel-1.png)

Figure 1: Shedding profile, relative to prodrome onset. The dashed line
marks rash onset.

</div>

### Parameters

The defaults of `measles_ww_params()`, with their sources, following the
format of the parameter table in the measles package
(`measles::measles_parameters()`). The table is stored in
[data/parameters.csv](data/parameters.csv), and the tests check that it
matches the defaults in the code. Values marked “Assumption” have no
direct source; the notes say what informed them and how they differ from
the measles package defaults.

<div id="tbl-params">

Table 1

<div class="cell-output-display">

| Parameter | Default | Literature range | Source | Notes |
|:---|:---|:---|:---|:---|
| Population size (`n`) | `1e4` agents |  | Scenario-specific input | Age distribution of Salt Lake County, UT (2025 Census estimates, via multigroup.vaccine), rescaled to n. |
| Initial infections (`n_seed`) | `1` agents |  | Assumption: a single introduction |  |
| R0 (`R0`) | `7` secondary cases | 12–18 | Assumption [\[1\]](https://doi.org/10.1016/S1473-3099(17)30307-9) | Below the 12-18 range of Guerra et al. 2017 (Lancet Infect Dis 17(12):e420-e428); the measles package uses 15. Not a model argument: it enters through the calibrated contact matrix. |
| Vaccination rate (`prop_vaccinated`) | `0` proportion |  | Scenario-specific input | The analysis runs 0% and 85% coverage. |
| Vax efficacy (`vax_efficacy`) | `0.97` per vaccinated agent |  | Utah DHHS Measles Disease Plan (updated March 26, 2026) (“~97%”) [\[1\]](https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf) | Same as the measles package. |
| Transmission rate (`transmission_rate`) | `0.9` per contact |  | Assumption: highly transmissible; Utah DHHS Measles Disease Plan (updated March 26, 2026): “90% of susceptible contacts will develop disease” [\[1\]](https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf) | Fixed high, then the contact matrix is calibrated to R0. Same as the measles package. |
| Incubation period (`incubation_period`) | `10` days | 8–12 | Utah DHHS Measles Disease Plan (updated March 26, 2026): exposure to prodrome averages 8-12 days [\[1\]](https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf) | The measles package uses 12. |
| Prodromal period (`prodromal_period`) | `4` days | 2–4 | Utah DHHS Measles Disease Plan (updated March 26, 2026): prodrome lasts 2-4 days (range 2-8); contagious 4 days before rash onset [\[1\]](https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf) | Same as the measles package. |
| Rash period (`rash_period`) | `4` days | 2–4 | Utah DHHS Measles Disease Plan (updated March 26, 2026): contagious to 4 days after rash onset [\[1\]](https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf) | The measles package uses 3 (infectivity is minimal after day 2 of rash). |
| Rash reduction contact rate (`rash_reduction_contact_rate`) | `0.5` proportion |  | Assumption | The measles package uses 1 (full reduction). |
| Probability of hospitalization (`prob_hospitalization`) | `0.2` per case | 0.08–0.185 | Assumption [\[1\]](https://doi.org/10.1056/EVIDpha2600141) [\[2\]](https://doi.org/10.15585/mmwr.mm7520a1) | Converted to the daily rate that ModelMeaslesMixing() takes. Observed: 8% overall in Utah (Jones et al. 2026, NEJM Evid 5(8)); 18.5% in West Texas (Wang et al. 2026, MMWR 75:252-257). Does not affect shedding or triggers. |
| Hospitalization period (`hospitalization_period`) | `7` days |  | Assumption [\[1\]](https://doi.org/10.1056/EVIDpha2600141) | Observed stays are shorter: a mean of 2.1 nights in Utah (Jones et al. 2026). |
| Simulation length (`ndays`) | `365` days |  | Scenario-specific input |  |
| Peak shedding (`ww_t_peak`) | `4` days after prodrome onset |  | Corchis-Scott et al. 2026, Emerg Infect Dis 32(10): shedding during the prodrome and around rash onset [\[1\]](https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article) |  |
| End of shedding (`ww_t_end`) | `18` days after prodrome onset |  | Corchis-Scott et al. 2026: urine RNA up to 14 days; CDC: RT-PCR positive up to 10+ days after rash [\[1\]](https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article) |  |
| Peak-to-tail ratio (`ww_log10_vpeak`) | `3` log10 |  | Assumption | Only sets the shape; the profile is normalized to sum to one. |
| Copies per infection, mean (`ww_log10_g_mean`) | `10` log10 genome copies |  | STOTEN 2024: about 10^8.7 copies per infection per day over about 2 weeks; Corchis-Scott et al. 2026: median 1.4 x 10^11 copies per person per day [\[1\]](https://www.sciencedirect.com/science/article/abs/pii/S0048969724082998) [\[2\]](https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article) |  |
| Copies per infection, SD (`ww_log10_g_sd`) | `1` log10 genome copies |  | Corchis-Scott et al. 2026: shedding spans four orders of magnitude [\[1\]](https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article) |  |
| Wastewater per person (`ww_alpha`) | `378500` mL per person per day |  | CDC wwinference model [\[1\]](https://github.com/CDCgov/ww-inference-model) |  |
| Measurement SD (`ww_sigma`) | `0.5` log scale |  | Assumption |  |
| Limit of detection (`ww_lod`) | `1` copies/mL |  | STOTEN 2024: process LOD [\[1\]](https://www.sciencedirect.com/science/article/abs/pii/S0048969724082998) |  |
| Sampling days (`ww_sampling_days`) | `c(1, 3, 5)` days of the week |  | Corchis-Scott et al. 2026 (Ontario) [\[1\]](https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article) | Utah samples measles twice a week (Deseret News, March 21, 2026: https://www.deseret.com/lifestyle/2026/03/21/measles-outbreak-utah-public-health-lab-wastewater-surveillance/). Sensitivity: 1, 2, 3, and 5 times a week. |
| Lab turnaround (`ww_lab_delay`) | `3` days | 2–10 | Assumption, informed by: Utah pilot (2020), transport within 24 h and a 48-h lab turnaround; Provo, UT (2025), 4 days from collection to DHHS notification by WastewaterSCAN; NWSS, a median of 2 days (IQR 1-4) from collection to testing; New Mexico (2025), results within 10 days [\[1\]](https://www.acwa-us.org/wp-content/uploads/2020/08/COVID-Surv._Utah_wastewater_virus_pilot_Aug.pdf) [\[2\]](https://dhhs.utah.gov/featured-news/measles-has-been-detected-in-provo-area-wastewater/) [\[3\]](https://pmc.ncbi.nlm.nih.gov/articles/PMC11103741/) [\[4\]](https://pmc.ncbi.nlm.nih.gov/articles/PMC13150633/) | The Utah lab figure is from the COVID-19 pilot; measles assays may take longer. Sensitivity: 0, 1, 3, 5, 7, and 10 days. |
| Case detection probability (`case_report_prob`) | `0.7` per case |  | Assumption; Corchis-Scott et al. 2026 estimated 5-10x under-reporting in Ontario [\[1\]](https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article) | Sensitivity: 70%, 50%, 30%, and 10%. |
| Case reporting delay (`case_report_delay`) | `3` days |  | Assumption | The measles package assumes 2 days from rash onset to detection. |
| Case trigger (`trigger_cases`) | `5` reported cases |  | Assumption |  |
| WW trigger (`trigger_ww`) | `1` copies/mL |  | Calibrated (see Calibrating the triggers) | The default is a placeholder; the analysis replaces it with the calibrated threshold. |
| School contact reduction (`school_contact_reduction`) | `0.8` proportion |  | Assumption |  |

</div>

</div>

### Calibrating the triggers

The clinical system sees a case at rash onset, but only after the
reporting delay. WW sees it from the prodrome, but only after the lab
turnaround. The WW threshold is calibrated on the baseline runs so that,
at the median, both triggers fire after the same number of cases:

- **Cases**: cumulative rash onsets on the rash day of the report that
  crosses the case threshold.
- **WW**: cumulative prodrome onsets on the collection day of the first
  sample above the WW threshold.

Both medians are taken over the same runs, those where the case trigger
fires. If WW never crosses the threshold in one of those runs, it counts
as WW never seeing the outbreak, so a threshold cannot look well
calibrated just by missing the outbreaks it would detect late. The
calibration uses the reference case detection probability (70%).

Because runs with the same seed share their trajectory until the
trigger, the scenarios are paired, and the lead time is the case trigger
day minus the WW trigger day in the same outbreak.

## Results

Each scenario was simulated 400 times (`analysis/run-scenarios.R`) for
two levels of vaccination coverage: none (R0 = 7 in a fully susceptible
city) and 85% (effective R ≈ 1.2). Each simulation starts from a single
introduction. Summaries only include runs in which the baseline produced
an outbreak (50 or more infections), and the “Outbreaks” column gives
their share. “Infected at trigger” and “Reported at trigger” are medians
on the trigger day.

<div id="tbl-summary">

Table 2

<div class="cell-output-display">

| Coverage | Outbreaks (%) | Scenario | WW threshold (copies/mL) | Trigger day | Infected at trigger | Reported at trigger | Infections, median \[IQR\] | Averted vs baseline, mean | Lead of WW over cases (days) |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 0% | 96 | baseline |  |  |  |  | 9,976 \[9,972, 9,979\] |  |  |
| 0% | 96 | wastewater | 4.5 | 20 | 57 | 3 | 9,952 \[9,948, 9,957\] | 23 | 3 |
| 0% | 96 | cases |  | 24 | 110 | 5 | 9,952 \[9,948, 9,956\] | 24 |  |
| 85% | 40 | baseline |  |  |  |  | 891 \[836, 938\] |  |  |
| 85% | 40 | wastewater | 11.2 | 36 | 13 | 4 | 184 \[31, 377\] | 656 | 5 |
| 85% | 40 | cases |  | 42 | 16 | 5 | 230 \[52, 369\] | 639 |  |

</div>

</div>

``` r
curves <- rbindlist(lapply(res, `[[`, "curves"))
curves[, coverage := paste0("Vaccination coverage: ", 100 * coverage, "%")]
ggplot(curves, aes(day, new_inf, colour = scenario)) +
  geom_line() +
  facet_wrap(~coverage, scales = "free_y") +
  labs(x = "Day", y = "New infections (mean)", colour = "Scenario") +
  coord_cartesian(xlim = c(0, 200)) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

<div id="fig-curves">

![](README_files/figure-commonmark/fig-curves-1.png)

Figure 2: Mean daily new infections by scenario.

</div>

### Takeaways

- **Without vaccination**, WW triggers a median of 3 days before the
  case-based trigger. At that point about 57 people are infected and 3
  cases have been reported. The case trigger fires at about 110 infected
  and 5 reported. Both thresholds respond to the same number of cases
  (about 8); WW simply sees them earlier. However, the school-age
  intervention only lowers R from 7 to about 5.4, so neither trigger
  stops a city-wide outbreak.
- **With 85% coverage** (effective R ≈ 1.2), the same intervention
  brings R just below one (≈ 0.95) and averts most infections. But the
  outbreak grows slowly, and at 70% case detection the WW lead is about
  5 days, so the differences between the WW and case-based scenarios are
  within simulation noise. WW gains a clear advantage only when clinical
  detection is lower (see [Case detection
  probability](#case-detection-probability)).
- The WW lead comes from shedding during the prodrome (about 4 days
  before rash) plus the case-reporting delay, minus the lab turnaround.
  Each extra day of turnaround removes about one day of lead. With a
  7-day turnaround, the median lead is -1 days without vaccination and 1
  days with 85% coverage (table below).
- Without vaccination, sampling once a week instead of 3 times a week
  cuts the lead from 3 to 0 days, about what the extra days between
  samples would predict. With 85% coverage it costs far more, likely
  because the slowly growing signal stays near the threshold for weeks,
  so each missed sampling day is a missed chance to cross it. The WW
  trigger moves from day 36 to day 53, and the lead from 5 to -3 days.
  Sampling twice a week, as in Utah, gives a median lead, in days, of 1
  without vaccination and 1 with 85% coverage (see [Sampling
  frequency](#sampling-frequency)).
- Earlier detection pays off most when the outbreak grows fast *and* the
  triggered intervention is strong enough to bring R below one. With R0
  = 7 and no vaccination, that would take broader measures than reducing
  school-age contacts.

### Lab turnaround time

The value of earlier detection shrinks as the lab turnaround time grows.
The default of 3 days is an assumption. It matches what Utah reports: in
the 2020 pilot, samples reached the Utah Public Health Laboratory within
24 hours and results took 48 hours, and in the July 2025 Provo
detection, DHHS was notified 4 days after collection. Other programs are
slower (results within 10 days in New Mexico), and measles assays may
take longer than the COVID-19 assays of the Utah pilot, so the range
goes up to 10 days. The WW threshold is fixed at its calibrated value:

<div id="tbl-lab">

Table 3

<div class="cell-output-display">

| Coverage | Lab turnaround (days) | Trigger day | Lead of WW over cases (days) | Infections averted, mean |
|:---|---:|:---|:---|:---|
| 0% | 0 | 17 | 6 | 23 |
| 0% | 1 | 18 | 5 | 23 |
| 0% | 3 | 20 | 3 | 23 |
| 0% | 5 | 22 | 1 | 24 |
| 0% | 7 | 24 | -1 | 23 |
| 0% | 10 | 27 | -4 | 23 |
| 85% | 0 | 33 | 8 | 654 |
| 85% | 1 | 34 | 7 | 648 |
| 85% | 3 | 36 | 5 | 656 |
| 85% | 5 | 38 | 3 | 643 |
| 85% | 7 | 40 | 1 | 631 |
| 85% | 10 | 43 | -2 | 624 |

</div>

</div>

### Sampling frequency

WW is sampled 3 times a week by default, as in Ontario; Utah samples
measles twice a week. Less frequent sampling delays the first sample
above the threshold and gives fewer chances to catch it. The WW
measurement error is drawn every day, so the schedules share their
measurements and differ only in which days are sampled. The threshold is
fixed at its calibrated value, and the lead is measured against the
case-based trigger in the same runs:

<div id="tbl-sampling">

Table 4

<div class="cell-output-display">

| Coverage | Samples per week | Triggered (%) | Trigger day | Infected at trigger | Lead of WW over cases (days) | Infections averted, mean |
|:---|---:|:---|:---|:---|:---|:---|
| 0% | 1 | 100 | 25 | 123 | 0 | 23 |
| 0% | 2 | 100 | 21 | 78 | 1 | 24 |
| 0% | 3 | 100 | 20 | 57 | 3 | 23 |
| 0% | 5 | 100 | 20 | 52 | 3 | 23 |
| 85% | 1 | 100 | 53 | 20 | -3 | 627 |
| 85% | 2 | 100 | 40 | 15 | 1 | 638 |
| 85% | 3 | 100 | 36 | 13 | 5 | 656 |
| 85% | 5 | 100 | 36 | 11 | 6 | 666 |

</div>

</div>

### Case detection probability

Measles cases are hard to detect clinically: many are never seen by a
clinician, tested, or reported. Here each case is detected during its
rash period with a given probability, and the WW threshold stays at the
value calibrated for 70% detection. WW does not depend on clinical
detection, so only the case-based trigger changes. Case detection uses
common random numbers (one draw per case), so the cases detected at a
lower probability are a subset of those detected at a higher one, and
the lead is paired with the same WW runs. “Triggered” is the share of
outbreaks in which the case-based trigger fired; the lead is computed
over outbreaks in which both triggers fired.

<div id="tbl-report">

Table 5

<div class="cell-output-display">

| Coverage | Case detection (%) | Triggered (%) | Trigger day | Infected at trigger | Lead of WW over cases (days) | Infections averted, mean |
|:---|:---|:---|:---|:---|:---|:---|
| 0% | 70 | 100 | 24 | 110 | 3 | 24 |
| 0% | 50 | 100 | 26 | 160 | 4 | 23 |
| 0% | 30 | 100 | 28 | 269 | 7 | 23 |
| 0% | 10 | 100 | 32 | 782 | 11 | 22 |
| 85% | 70 | 100 | 42 | 16 | 5 | 639 |
| 85% | 50 | 100 | 49 | 20 | 9 | 606 |
| 85% | 30 | 100 | 64 | 31 | 21 | 568 |
| 85% | 10 | 100 | 97 | 78 | 50 | 453 |

</div>

</div>

- Lower clinical detection delays the case-based trigger, so the WW lead
  grows: without vaccination it goes from 3 days at 70% detection to 11
  days at 10%.
- With 85% coverage the outbreak grows slowly, so missed cases cost more
  time: the lead goes from 5 days at 70% detection to 50 days at 10%,
  and the infections averted by the case-based trigger fall from 639 to
  453 (WW: 656).
- Without vaccination the averted infections barely change, because the
  intervention cannot stop the outbreak regardless of when it starts.

The grid of detection probabilities is set by `report_probs` (and the
reference value by `report_prob_ref`) in `analysis/run-scenarios.R`.

## Reproducing

The analysis needs epiworldR 0.19.0, measles 0.4.1, and
multigroup.vaccine 0.1.1, which are installed from GitHub at pinned
commits. The easiest way to get them is the dev container
(`.devcontainer/Containerfile`), which also pins the R version and the
CRAN snapshot:

- **In the browser**: [![Open in GitHub
  Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/EpiForeSITE/wastewater-econ-abm?devcontainer_path=.devcontainer/web/devcontainer.json)
  (uses `.devcontainer/web/devcontainer.json`).
- **Locally**: open the project in VS Code with the Dev Containers
  extension (uses `.devcontainer/devcontainer.json`), or build the image
  directly with
  `docker build -f .devcontainer/Containerfile .devcontainer`.

The image is built and tested by the `Dev container` GitHub Actions
workflow, and published as `ghcr.io/epiforesite/wastewater-econ-abm-dev`
from `main`. Inside the container:

``` r
source("analysis/run-scenarios.R")   # writes analysis/results.rds
quarto::quarto_render("README.qmd")
tinytest::run_test_dir("tests")
```

All parameters are set in `measles_ww_params()` (R/simulate.R) and can
be changed, e.g., `measles_ww_params(R0 = 12, ww_lab_delay = 5)`.

## References

- CDC wwinference model: <https://github.com/CDCgov/ww-inference-model>
- Corchis-Scott R. et al. Wastewater surveillance to track resurgent
  measles outbreak, Ontario, Canada, 2025. *Emerg Infect Dis*
  2026;32(10). <https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article>
- Temporal, spatial, and methodological considerations in evaluating the
  viability of measles wastewater surveillance. *Sci Total
  Environ* 2024.
  <https://www.sciencedirect.com/science/article/abs/pii/S0048969724082998>
- Measles RNA detection in wastewater solids: longitudinal monitoring at
  a national scale.
  <https://pubs.rsc.org/ew/article/12/7/2074/1223010/Measles-RNA-detection-in-wastewater-solids>
- Genotype-specific detection of measles virus using wastewater
  surveillance, Wisconsin, 2026. *MMWR* 75(27).
  <https://www.cdc.gov/mmwr/volumes/75/wr/mm7527a1.htm>
- Wastewater surveillance for measles virus during a measles outbreak,
  Colorado, 2025.
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12806976/>
- Utah DHHS. Measles Disease Plan (updated March 26, 2026).
  <https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf>
- Utah DHHS. Measles has been detected in Provo area wastewater (July
  2025).
  <https://dhhs.utah.gov/featured-news/measles-has-been-detected-in-provo-area-wastewater/>
- Why Utah Public Health Lab uses sewage to track measles. *Deseret
  News*, March 21, 2026.
  <https://www.deseret.com/lifestyle/2026/03/21/measles-outbreak-utah-public-health-lab-wastewater-surveillance/>
- Use of wastewater surveillance in the coordinated response to COVID-19
  in Utah (pilot report, August 2020).
  <https://www.acwa-us.org/wp-content/uploads/2020/08/COVID-Surv._Utah_wastewater_virus_pilot_Aug.pdf>
- CDC NWSS sampling sites, time from collection to testing.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11103741/>
- Measles wild-type virus detection through wastewater surveillance in
  Sandoval County, New Mexico.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC13150633/>
- Guerra F. M. et al. The basic reproduction number (R0) of measles: a
  systematic review. *Lancet Infect Dis* 2017;17(12):e420-e428.
  <https://doi.org/10.1016/S1473-3099(17)30307-9>
- Litvinova M. et al. Epistorm-Mix. medRxiv 2025.11.20.25340662.
- Klinkenberg D., Nishiura H. The correlation between infectivity and
  incubation period of measles. *J Theor Biol* 2011.
