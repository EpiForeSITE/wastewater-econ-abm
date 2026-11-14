# Project guidance

- This is a measles wastewater project. The main `build_model()` and scenario
  runners must use the pinned `epiworldR` `v1-layers` implementation, with
  school and community contact layers sharing agents. A separate example alone
  is not a migration.
- The current application uses Salt Lake County population and school
  vaccination data; keep these choices in application adapters, with documented
  source versions, geographic and enrollment assumptions. Coverage is not the
  same as vaccine effectiveness.
- Keep `README.qmd`, its rendered `README.md`, scripts, saved results, and
  figures consistent. Clearly label older simulation results if they have not
  been rerun; do not present them as outputs of the new model.
- Document natural-history parameters and assumptions using roxygen2 comments,
  distinguish source evidence for each disease stage, and test the actual main
  scenario path.
- Keep the analysis simple: use one short R helper script and run scenarios
  directly from the Quarto document with `epiworldR::run_multiple()`. Keep
  tests in their own folder; avoid wrapper chains and unnecessary runners.
- For this work, open pull requests against `measles-wastewater-scenarios`
  unless the user specifies otherwise.

## R implementation conventions

- Prefer the native `|>` pipe for sequential model-builder operations. Load
  epiworldR once at the entry point rather than repeating `epiworldR::`
  throughout builder code.
- Apply the 80-column limit and multiline-call formatting across all project
  files, including `stopifnot()` and `c()`: break after `(` and before `)`, name
  arguments where meaningful, and align `=` signs.
- Separate logical code sections with blank lines and brief `#` comments
  explaining their purpose.
- Keep reusable functions independent of the application location. Utah and Salt
  Lake County belong in application adapters and scripts; accept population,
  school and contact inputs for other locations.
- Modeling and population functions consume the caller's RNG stream. Do not set
  or save/restore R seeds inside them. Initialize seeds in execution scripts,
  including explicit paired-replicate boundaries.
- Use Git history for earlier reports and parameter versions. Do not add
  historical copies to the working tree.

- When asked to address review comments, complete the requested changes, run
  relevant checks, commit and push, and reply directly to each comment unless
  the user says otherwise.
