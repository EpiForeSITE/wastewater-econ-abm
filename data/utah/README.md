# Utah school data

`school_district_schoolID_data_Utah_2025.csv` is copied from
[`multigroup-vaccine/data-raw`](https://github.com/EpiForeSITE/multigroup-vaccine/blob/ddd9201f28b9908c033542f6a2445b84619af1c7/data-raw/school_district_schoolID_data_Utah_2025.csv),
commit `ddd9201f28b9908c033542f6a2445b84619af1c7`.

The file contains school identifiers, district, county, level, `popfrac`,
`vaxrate`, and `priorrate`. In Salt Lake County, `popfrac` sums to about one
within each district, not within the county. Do not interpret it as a county
share or a school enrollment count. The source model uses district census
estimates to allocate populations and interprets `vaxrate` as MMR2 coverage.
School IDs may occur at multiple levels; each ID and level pair is a record.
The exact source and year of the underlying vaccination and enrollment
measurements still need verification before scientific use.

`acs_Utah_2023_school_districts.csv` is the corresponding district ACS data
from the same [upstream commit](https://github.com/EpiForeSITE/multigroup-vaccine/blob/ddd9201f28b9908c033542f6a2445b84619af1c7/inst/extdata/acs_Utah_2023_school_districts.csv).
The model uses its ages 5--17 estimates to weight district-level school
fractions. These are residential population estimates, not school enrollment;
district boundaries and attendance do not necessarily coincide.
