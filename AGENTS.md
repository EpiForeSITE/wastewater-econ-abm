# Project guidance

- This is a measles wastewater project. The main `build_model()` and scenario runners must use the pinned `epiworldR` `v1-layers` implementation, with school and community contact layers sharing agents. A separate example alone is not a migration.
- Use Salt Lake County population and school vaccination data with documented source versions, geographic and enrollment assumptions. Coverage is not the same as vaccine effectiveness.
- Keep `README.qmd`, its rendered `README.md`, scripts, saved results, and figures consistent. Clearly label older simulation results if they have not been rerun; do not present them as outputs of the new model.
- Document natural-history parameters and assumptions using roxygen2 comments, distinguish source evidence for each disease stage, and test the actual main scenario path.
- For this work, open pull requests against `measles-wastewater-scenarios` unless the user specifies otherwise.
