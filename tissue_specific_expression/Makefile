#!/usr/bin/make

include ../mk/variables.mk

## eventually, this file should be downloadable from somewhere like
## dryad; in the meantime, it will be made available from one of my
## servers once I publish the abstract
combined_read_counts:
	if [ -e ../tissue_expression/combined_read_counts ]; then \
		ln -s ../tissue_expression/combined_read_counts $@; \
    else \
		wget -O $@ "http://some-place-to-be-determined/combined_read_counts"; \
    fi;

chosen_samples:
	if [ -e ../tissue_expression/chosen_samples ]; then \
		ln -s ../tissue_expression/chosen_samples $@; \
    else \
		wget -O $@ "http://some-place-to-be-determined/chosen_samples"; \
    fi;

GTEx_Data_V4_Annotations_SampleAttributesDS.txt:
	if [ -e ../tissue_expression/$@ ]; then \
		ln -s ../tissue_expression/$@ $@; \
    else \
		wget -O $@ "http://www.gtexportal.org/static/datasets/gtex_analysis_v4/annotations/$@" || \
			rm $@; \
	fi;



categorized_samples: R/categorize_samples.R chosen_samples \
	GTEx_Data_V4_Annotations_SampleAttributesDS.txt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

tissue_specific_markers: R/tissue_specific_expression.R categorized_samples \
	combined_read_counts
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@
