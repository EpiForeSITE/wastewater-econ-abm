# Render one report per scenario in scenarios.yml, then compare them.
# Results go to results/<scenario>/ and reports to scenarios/<scenario>.md.
# WASTEWATER_NSIMS and WASTEWATER_CORES set replicates and cores (1000, 5).

SCENARIOS := baseline wastewater wastewater-fast

.PHONY: all $(SCENARIOS) compare test clean
.NOTPARALLEL:

all: $(SCENARIOS) compare

$(SCENARIOS): %: scenarios/%.md

# Render a copy named after the scenario, so its figures do not collide.
scenarios/%.md: report.qmd model.R scenarios.yml data/parameters.csv
	cp report.qmd $*.qmd
	quarto render $*.qmd -P scenario:$* || (rm -f $*.qmd; exit 1)
	rm -f $*.qmd
	mkdir -p scenarios
	rm -rf scenarios/$*_files
	mv $*.md scenarios/$*.md
	if [ -d $*_files ]; then mv $*_files scenarios/$*_files; fi

# Paired comparison across scenarios.
compare: scenarios/README.md

scenarios/README.md: scenarios/README.qmd \
	$(foreach s,$(SCENARIOS),scenarios/$(s).md)
	quarto render scenarios/README.qmd

test:
	Rscript -e 'r <- tinytest::run_test_dir("tests"); print(r)' \
		-e 'if (!all(as.logical(r))) quit(status = 1)'

clean:
	rm -rf $(foreach s,$(SCENARIOS),scenarios/$(s).md scenarios/$(s)_files \
		results/$(s))
