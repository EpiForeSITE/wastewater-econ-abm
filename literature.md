# Literature notes

What we found in the literature about measles shedding, wastewater
surveillance and outbreak response, and how it shapes the model's
parameters ([data/parameters.csv](data/parameters.csv)) and scenarios
([scenarios.yml](scenarios.yml)). Numbered references are listed at the
end. Some sources were behind paywalls or blocked automated access. Where a
claim rests on an abstract or a secondary summary, the text says so.

## Natural history and shedding

- **Incubation.** About 10–12 days from infection to fever and about 14 to
  rash. The virus first replicates in the respiratory tract and lymphoid
  tissue, and infectious virus appears in the blood 4–7 days after
  infection [1, 2]. In experimentally infected rhesus macaques, viremia
  starts 7–10 days after infection and the rash appears on days 10–11 [3].
- **Infectious period.** Contagiousness runs from about 4 days before to
  4 days after rash onset. Utah [4] and Colorado [11] both use this window,
  and the model uses it: the prodrome lasts 4 days and the infectious rash
  4 days.
- **When shedding starts.** People shed measles RNA during the prodrome,
  before rash onset, mainly in urine [5]. We found no data on urine or
  fecal shedding before the prodrome, so shedding halfway through the
  latent period (about 9–10 days before rash) has no support. Shedding
  plausibly starts at prodrome onset or a day or two earlier. The model
  starts it at
  prodrome onset (`ww_t_start = 0`) and can start it up to `ww_lab_delay`
  days earlier.
- **How long and how much.** Urine RNA lasts about 1–2 weeks after rash,
  and RNA can persist at low levels for weeks to months [3, 5]. Corchis-Scott
  et al. estimated a median of 1.4 × 10^11 genome copies per person per day
  from hospital effluent, spanning four orders of magnitude [5]. The model
  sheds for 18 days from prodrome onset, peaking at rash onset, with a
  lognormal total per infection.
- **Vaccine strain.** Vaccinees shed vaccine-strain RNA in urine for up to
  14 days [5], which genotyping can distinguish from wild type. Such
  shedding is one source of false positives in assays that do not
  genotype; the model uses 0.5% per sample.

## Wastewater surveillance in practice

**Sampling and turnaround.** Utah tests for measles twice a week [6]. In
2025, a Provo sample reached Utah DHHS 4 days after collection [7]. The
national system reports a median of 2 days (interquartile range 1–4) from
collection to testing [8], and New Mexico reported results within 10 days
[9]. A 2020 Utah pilot shipped samples within 24 hours and had a 48-hour
lab turnaround [10]. The `wastewater` scenario uses 4 days and
`wastewater-fast` 2.

**Lead over case reports.** Field reports disagree, and measure different
things:

| Place | Wastewater lead | Compared with | Source |
|:---|:---|:---|:---|
| Mesa County, Colorado, 2025 | Samples 1–3 days before the index case's symptom onset | Symptom onset; tied to an unrecognized family cluster | [11] |
| Sandoval County, New Mexico, 2025 | 5 days | Clinical case detection | [9] |
| Texas, 2025 | About 10–14 days | Official outbreak announcement | [12] |
| Oregon, 2024 | 10 weeks | First confirmed case | [13] |
| Ontario, 2025 | None; wastewater lagged weekly cases by about a week | Weekly cases | [5] |
| Western Switzerland, 2024 | None; tracked rash onset | Rash onset | [14] |

Leads of weeks (Oregon, Texas) cannot come from shedding before rash, which
lasts at most about a week. They reflect transmission nobody detected: most
cases go unreported. The model reproduces this. Its lead comes from
under-reporting, not from early shedding.

**What health departments do after a detection.** Responses have been
alerts and vaccine access, not mass testing [15]:

- Southern Nevada (August 2025): a detection with no reported cases. The
  health district urged families to check vaccination status and extended
  clinic hours.
- St. Mary's County, Maryland (June 2025): added measles to wastewater
  monitoring and scheduled extra vaccine clinics.
- Hawaiʻi (August 2025) and Georgia (October 2025): health advisories asking
  clinicians to test patients with compatible symptoms and to promote
  vaccination.
- Mecklenburg County, North Carolina (January 2026): walk-in MMR at all
  county immunization clinics.
- Utah DHHS issued a public notice after the Provo detection [7].

Officials stress that a single detection does not confirm cases. This
supports the model's design: a positive result starts the same alert as a
reported case, and the alert includes a vaccination campaign.

## Case-based surveillance and under-reporting

- **Reporting.** Genomic and transmission models of the 2025–26 southwest
  Utah outbreak estimate that about 25% of cases were reported [16]. In
  Ontario, wastewater implied 5–10 times more cases than were reported [5].
  The model reports 25% of infections outside alerts.
- **Hospitalization.** 8% of reported cases in Utah were hospitalized [17].
  Times 25% reporting, that is 2% of infections. West Texas reported 18.5%
  [18].
- **During an alert.** We found no estimate of how much reporting rises
  once an outbreak is known. Health alerts to clinicians, testing and
  exposure notification all raise it, and the model assumes 50% during an
  alert (`alert_reporting_fraction`).

## Outbreak response

- **Alert and duration.** An outbreak is declared over after two incubation
  periods (42 days) with no new cases. The model keeps an alert on for
  42 days after the last trigger.
- **School exclusion.** Exposed students without evidence of immunity who
  do not receive post-exposure prophylaxis are excluded for 21 days after
  their last exposure [19, 20]. The model excludes a school's unvaccinated
  students for 42 days after the school's last reported case, which is
  longer. Exposure dates are not tracked.
- **Post-exposure prophylaxis.** MMR within 72 hours of exposure, or immune
  globulin within 6 days [21]. In the 2013 New York City outbreak, MMR
  prophylaxis was 83.4% effective (95% CI 34.4–95.8%) and immune globulin
  100% (56.2–99.8%) [21]; the intervals are wide, with few events.

## Vaccination campaigns

Each outbreak response reviewed here included an MMR campaign: clinics,
school sessions, outreach through community leaders, and early doses for
infants aged 6–11 months.

### How much uptake rises

| Outbreak | Response | Uptake | Source |
|:---|:---|:---|:---|
| Minnesota, 2017 (75 cases) | Accelerated MMR recommendation for Somali Minnesotan children and outbreak counties; culturally targeted outreach | Statewide weekly doses rose from about 2,700 to 9,964 (3.7 times); more than 51,000 doses above the expected baseline | [19, 20] |
| Clark County, Washington, 2019 (71 cases) | Clinics and messaging | MMR orders in January up nearly 500% (530 to 3,150); Kaiser clinics up 417% | [22] |
| New York City, 2018–19 (649 cases) | Mandatory vaccination order in four ZIP codes (April 2019), outreach | Williamsburg coverage of children aged 12–59 months rose from 79.5% to 91.1% in about 11 months | [23, 24] |
| Swansea, Wales, 2013 (447 cases) | Weekend hospital clinics, then school sessions | 60% of 9,000 unvaccinated 10–18-year-olds received a dose by the end of May, about 8 weeks; more than 77,000 doses overall | [25] |
| West Texas, 2025 (762 cases) | Clinics, early infant doses | Region 1 registry doses up 60% in January–March (5,104 to 8,189); statewide up 16% | [18, 26] |
| New Mexico, 2025 | Outbreak response | MMR doses up 55% (January–September, year on year) | [9] |
| Chicago shelter, 2024 (57 cases) | Mass vaccination at the shelter within 1 day of confirming the index case | 93% coverage within 4 days | [27] |
| Southwest Utah, 2025–26 | Clinics, early infant doses | A "slight uptick", mostly hesitant parents choosing one dose | [28] |

### The model's campaign rate

Increases in doses are relative to a baseline, so they do not give the
share of unvaccinated people reached. Two responses do:

- **New York City:** the unvaccinated share fell from 20.5% to 8.9% in
  about 335 days, a daily rate of ln(20.5 / 8.9) / 335 ≈ 0.0025. This
  included a vaccination mandate, and new birth cohorts dilute it.
- **Swansea:** 60% of the unvaccinated reached in about 56 days, a daily
  rate of ln(1 / 0.4) / 56 ≈ 0.016, through school sessions.

The model uses `campaign_rate = 0.005`: each day of an alert, every agent
not yet vaccinated (and not detected as a case) gets a dose with
probability 0.005. That is about 19% of the unvaccinated over one 42-day
alert, between New York City and Swansea. Southwest Utah's modest response
argues against the upper end. This is the parameter to vary in a
sensitivity analysis.

### When campaigns start

Campaigns start as soon as an outbreak is recognized: within a day of
confirming the index case in Chicago [27]. But recognition can be late:
Minnesota's first two cases were confirmed when the outbreak was already
in its third generation [20]. Chicago shows what
timing is worth: a model fitted in real time estimated that a week's delay
in mass vaccination would have raised the chance of 100 or more cases from
about 1% to 8% (15% without active case-finding), and that starting a week
earlier would have prevented it [27]. In a school outbreak model, responses
within 1–3 weeks prevented large outbreaks and responses within 6–7 weeks
still reduced them [29]. Recent wastewater detections led to vaccine
clinics before any case was reported [15]. The model starts the campaign
on the day the alert starts.

### How fast doses protect

MMR within 72 hours of exposure prevents about 83% of cases [21], so
vaccine immunity develops faster than infection from an exposure on the
same day. The model protects a susceptible agent from the day after the
dose, with the usual `vax_efficacy`; agents already infected are not
protected.

## Gaps

- No day-by-day profile of urine or fecal shedding for wild-type
  infection; the shedding kernel's shape is an assumption.
- No estimate of reporting during an alert.
- No campaign uptake data for Salt Lake County, or for any U.S. population
  that reports the share of unvaccinated people reached per day without a
  mandate.
- Little evidence on whether a wastewater detection, without cases, moves
  hesitant families to vaccinate.
- Not read in full: the WHO summary on environmental surveillance for
  measles, the national wastewater solids study (medRxiv
  10.1101/2025.07.18.25331801), and the Swiss study (abstract only).

## References

1. Review of measles pathogenesis and clinical course (read in an earlier
   search; bibliographic details not rechecked).
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC7149578/>
2. Griffin DE, et al. Measles virus, immune control, and persistence.
   *FEMS Microbiol Rev*. 2012 (PMID 22316382).
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC3319515/>
3. Association of persistent wild-type measles virus RNA with long-term
   humoral immunity in rhesus macaques. *JCI Insight*.
   <https://doi.org/10.1172/jci.insight.134992>
4. Utah DHHS. Measles disease plan (updated March 26, 2026).
   <https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf>
5. Corchis-Scott R, et al. Wastewater surveillance to track resurgent
   measles outbreak, Ontario, Canada, 2025. *Emerg Infect Dis*.
   2026;32(10). <https://wwwnc.cdc.gov/eid/article/32/10/26-0092_article>
6. Deseret News. Utah public health lab and wastewater surveillance of
   measles, March 21, 2026.
   <https://www.deseret.com/lifestyle/2026/03/21/measles-outbreak-utah-public-health-lab-wastewater-surveillance/>
7. Utah DHHS. Measles has been detected in Provo-area wastewater (2025).
   <https://dhhs.utah.gov/featured-news/measles-has-been-detected-in-provo-area-wastewater/>
8. National Wastewater Surveillance System turnaround times.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC11103741/>
9. Measles outbreak — New Mexico, 2025. *MMWR*. 2026;75(9).
   <https://www.cdc.gov/mmwr/volumes/75/wr/mm7509a1.htm>; and measles
   wild-type virus detection through wastewater surveillance in Sandoval
   County, New Mexico. *JAMA Netw Open*.
   <https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2848680>
10. Utah wastewater pilot (2020): transport within 24 hours and a 48-hour
    lab turnaround.
    <https://www.acwa-us.org/wp-content/uploads/2020/08/COVID-Surv._Utah_wastewater_virus_pilot_Aug.pdf>
11. Jensen et al. Notes from the field: wastewater surveillance for measles
    virus during a measles outbreak — Colorado, August 2025. *MMWR*.
    2026;75(2):20–22. <https://pmc.ncbi.nlm.nih.gov/articles/PMC12806976/>
12. Detection of measles in Texas wastewater. medRxiv, 2025.
    <https://www.medrxiv.org/content/10.1101/2025.04.08.25325475v1.full>
13. Notes from the field: retrospective analysis of wild-type measles virus
    in wastewater — Oregon, 2024. *MMWR*.
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12806952/>
14. Retrospective wastewater tracking of a measles outbreak in western
    Switzerland in winter 2024. *Environ Sci Technol Lett*.
    2025;12(6):689.
    <https://pubs.acs.org/doi/10.1021/acs.estlett.5c00159>
15. Health department releases on wastewater detections: Southern Nevada
    Health District (August 2025),
    <https://www.southernnevadahealthdistrict.org/news-release/measles-virus-detected-in-clark-county-wastewater/>;
    St. Mary's County, Maryland (June 2025),
    <https://smnewsnet.com/?p=548105>; Hawaiʻi DOH (August 2025),
    <https://health.hawaii.gov/news/?p=9775>; Georgia DPH (October 2025),
    <https://dph.georgia.gov/document/document/health-alert-measles-detection-wastewater-october-2025-0/download>;
    Mecklenburg County (January 2026),
    <https://cardinalpine.com/2026/01/13/measles-detected-in-charlotte-wastewater/>.
16. EpiForeSITE. Hidden burden of a measles outbreak revealed by genomic
    and transmission models.
    <https://github.com/EpiForeSITE/measles-hidden-burden>
17. Jones et al. *NEJM Evid*. 2026;5(8).
    <https://doi.org/10.1056/EVIDpha2600141>
18. Wang et al. Measles outbreak — West Texas, 2025. *MMWR*.
    2026;75:252–257. <https://doi.org/10.15585/mmwr.mm7520a1>; and
    <https://www.cdc.gov/mmwr/volumes/75/wr/mm7523a2.htm>
19. Hall V, et al. Measles outbreak — Minnesota, April–May 2017. *MMWR*.
    2017;66(27):713–717.
    <https://www.cdc.gov/mmwr/volumes/66/wr/mm6627a1.htm>
20. Banerjee E, et al. Containing a measles outbreak in Minnesota, 2017:
    methods and challenges. *Perspect Public Health*. 2020;140(3):162–171
    (abstract). <https://eld.bl.uk/catalog/100102431685_0x000001>
21. Arciuolo RJ, et al. Effectiveness of measles vaccination and immune
    globulin post-exposure prophylaxis in an outbreak setting — New York
    City, 2013. *Clin Infect Dis*. 2017;65(11):1843–1847.
    <https://doi.org/10.1093/cid/cix639>
22. The Columbian. MMR vaccines skyrocket locally, February 8, 2019.
    <https://columbian.com/news/2019/feb/08/mmr-vaccines-skyrocket-locally>;
    and NBC News, February 2019.
    <https://www.nbcnews.com/health/kids-health/measles-outbreak-sends-vaccine-demand-soaring-even-among-hesitant-n968806>
23. Zucker JR, et al. Consequences of undervaccination — measles outbreak,
    New York City, 2018–2019. *N Engl J Med*. 2020;382:1009–1017.
    <https://doi.org/10.1056/NEJMoa1911514>
24. Rosen JB. Measles in the Northeast (NYC DOHMH presentation, 2019).
    <https://www.maine.gov/dhhs/mecdc/infectious-disease/documents/2019-speaker-presentations/7b_MeaselsinNortheastRosen.pdf>
25. Public Health Wales. Measles surveillance and epidemiology.
    <https://phw.nhs.wales/topics/immunisation-and-vaccines/immunisation-surveillance/measles-surveillance-and-epidemiology>;
    and Senedd Health and Social Care Committee evidence, 2013.
    <https://business.senedd.wales/documents/s19074/Paper%201.html?CT=2>
26. Texas Medical Association. Vaccine uptake improves amid measles
    momentum (2025). <https://www.texmed.org/TexasMedicineDetail.aspx?id=65856>
27. Gressick K, et al. Measles outbreak associated with a migrant shelter —
    Chicago, Illinois, February–May 2024. *MMWR*. 2024;73(19).
    <https://www.cdc.gov/mmwr/volumes/73/wr/mm7319a1.htm>; and Real-time
    use of a dynamic model to measure the impact of public health
    interventions on measles outbreak size and duration — Chicago,
    Illinois, 2024. *MMWR*. 2024;73(19).
    <https://www.cdc.gov/mmwr/volumes/73/wr/mm7319a2.htm>
28. Salt Lake Tribune, April 22, 2026.
    <https://www.sltrib.com/news/politics/2026/04/22/utah-measles-outbreak-how-health/>;
    and KFF Health News.
    <https://kffhealthnews.org/public-health/utah-measles-outbreak-vaccines-preventable-diseases-doctors-strained-new-normal/>
29. Bonačić Marinović AA, et al. Effectiveness and timing of vaccination
    during school measles outbreak. *Emerg Infect Dis*. 2012;18(9):1405–1413.
    <https://wwwnc.cdc.gov/eid/article/18/9/11-1578>
