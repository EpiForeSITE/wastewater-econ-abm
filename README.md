# Measles wastewater surveillance in Salt Lake County

This repository simulates measles outbreaks in a synthetic Salt Lake County
population to estimate what a wastewater surveillance program adds to the
usual case-based response. It produces outcomes for an economic evaluation:
infections, hospital admissions and bed-days, school exclusions, wastewater
samples, public health alerts, and vaccination campaign doses. The model
is an
exploratory simulation, not a forecast.

## The model

The model is an agent-based model built with the `v1-layers` branch of
[epiworldR](https://github.com/UofUEpiBio/epiworldR) (pinned in the
devcontainer). It has 120,000 agents, a tenth of the county, with its age
structure and school MMR coverage.

- **Contacts.** Everyone mixes in a community layer by age group. Students
  also join one school of about 500, where every pair of schoolmates meets
  with the same daily probability. Contacts come from the Epistorm-Mix
  setting matrices; school contacts are moved from the community layer to the
  school layer, so a student's total contacts match the data.
- **Disease.** Infection, a latent period, an infectious prodrome, and an
  infectious rash. A quarter of cases are reported, as estimated for an
  outbreak in southwest Utah: hospitalized cases (2% of infections)
  and cases detected during rash are isolated and reported three days
  later.
- **Public health alert.** Any reported case starts a 42-day county-wide
  alert (extended by each new trigger). During it, half of cases are
  reported instead of a quarter, and a vaccination campaign gives each
  unvaccinated agent a 0.5% daily chance of an MMR dose (about a fifth of
  them over one alert).
- **School exclusions.** A reported student case excludes the school's
  unvaccinated students for 42 days after its last reported case.
- **Wastewater.** Infected agents shed into one county-wide catchment
  from prodrome onset (about four days before rash) for about two weeks.
  Samples are taken twice a week and return a result some days later. In
  the wastewater scenarios a positive result also starts the alert, so the
  program's value is starting the alert and its campaign earlier, net of
  false alarms.
- **Calibration.** The per-contact transmission probability is set so that
  an imported case causes seven secondary infections on average in a fully
  susceptible population.

```mermaid
flowchart LR
  subgraph contacts [Contact layers]
    C[Community<br/>everyone, by age]
    Sch[Schools<br/>students, about 500 each]
  end

  subgraph disease [Disease]
    S[Susceptible] --> L[Latent] --> P[Prodromal] --> R[Rash]
    R -->|"detected<br/>(25% reported)"| I[Isolated]
    R --> H[Hospitalized]
    I --> H
    R --> Rec[Recovered]
    I --> Rec
    H --> Rec
  end

  C -. infection .-> S
  Sch -. infection .-> S

  P --> Shed[Shedding into<br/>wastewater]
  R --> Shed
  Shed --> Sample[Sample<br/>twice a week] -->|lab delay| Pos{Positive?}
  Pos -->|"yes (wastewater<br/>scenarios)"| Alert["Public health alert, 42 days:<br/>50% of cases reported,<br/>vaccination campaign"]
  Alert -.-> I
  Alert -. MMR doses .-> S

  I -->|report delay| Alert
  I -->|report delay| Excl["School exclusion, 42 days:<br/>unvaccinated students<br/>leave the school layer"]
  Excl -.-> Sch
```

## Scenarios

Scenarios are defined in [scenarios.yml](scenarios.yml); each changes some
of the default parameters in [data/parameters.csv](data/parameters.csv).

| Scenario | Response |
|:---|:---|
| [baseline](scenarios/baseline.md) | Reported cases start alerts (with a vaccination campaign) and school exclusions |
| [wastewater](scenarios/wastewater.md) | Baseline, plus positive wastewater results start alerts: samples twice a week, results in 4 days |
| [wastewater-fast](scenarios/wastewater-fast.md) | As above, with results in 2 days |

Every scenario runs the same replicates: the same origin school, population,
vaccination and random draws, so outcomes can be compared replicate by
replicate. [scenarios/README.md](scenarios/README.md) compares them.
[literature.md](literature.md) summarizes the evidence behind the
shedding, surveillance and response assumptions.

## Workflow

1. Edit parameters in [data/parameters.csv](data/parameters.csv) or add a
   scenario to [scenarios.yml](scenarios.yml) (and to `SCENARIOS` in the
   [Makefile](Makefile)).
2. Run `make <scenario>` for one scenario, or `make` for all of them and the
   comparison. Each run copies [report.qmd](report.qmd) to
   `<scenario>.qmd`, renders it with the scenario as a parameter, and moves
   the report to `scenarios/<scenario>.md`.
3. Read the report in [scenarios/](scenarios), and use the tables in
   `results/<scenario>/`, described in [results/README.md](results/README.md).

`make test` runs the tests in [tests/](tests). `WASTEWATER_NSIMS` (default
1,000) and `WASTEWATER_CORES` (default 5) set the replicates and cores, e.g.
`WASTEWATER_NSIMS=10 make baseline` for a quick check. With 1,000
replicates on five cores, `make` takes about an hour and a half.

## Repository

| Path | Contents |
|:---|:---|
| [model.R](model.R) | Model builder, response event, and population helpers |
| [report.qmd](report.qmd) | Parameterized report that runs one scenario |
| [scenarios.yml](scenarios.yml) | Scenario definitions |
| [literature.md](literature.md) | Literature notes on shedding, surveillance and outbreak response |
| [data/](data) | Parameters, contact matrices, and Utah school data |
| [scenarios/](scenarios) | Rendered reports and the comparison |
| [results/](results) | Simulation outputs for the economic evaluation |
| [tests/](tests) | tinytest tests |
