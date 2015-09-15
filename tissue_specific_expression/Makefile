#!/usr/bin/make

include ../mk/variables.mk

VPATH=R

all: svm_classification tissue_specific_markers

## eventually, this file should be downloadable from somewhere like
## dryad; in the meantime, it will be made available from one of my
## servers once I publish the abstract
combined_read_counts:
	if [ -e ../tissue_expression/combined_read_counts ]; then \
		ln -s ../tissue_expression/combined_read_counts $@; \
    else \
		wget -O $@ "http://some-place-to-be-determined/combined_read_counts"; \
    fi;

combined_star_logs:
	if [ -e ../tissue_expression/$@ ]; then \
		ln -s ../tissue_expression/$@ $@; \
    else \
		wget -O $@ "http://some-place-to-be-determined/$@"; \
    fi;

combined_isoform_reads:
	if [ -e ../tissue_expression/$@ ]; then \
		ln -s ../tissue_expression/$@ $@; \
    else \
		wget -O $@ "http://some-place-to-be-determined/$@"; \
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

EG.name.txt:
	if [ -e ../tissue_expression/$@ ]; then \
		ln -s ../tissue_expression/$@ $@; \
    else \
		wget -O $@ "http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/$@" || \
			rm $@; \
	fi;

%_test_results: %_training.R svm_testing.R svm_classification
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< $(wordlist 3,$(words $^),$^) $@

%_test_results_reduced: %_training.R svm_testing_reduced.R svm_classification
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< $(wordlist 3,$(words $^),$^) $@

%_classification: %_training.R classification_svm.R interesting_gene_reads
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< $(wordlist 3,$(words $^),$^) $@

categorized_samples: R/categorize_samples.R chosen_samples \
	GTEx_Data_V4_Annotations_SampleAttributesDS.txt \
	EG.name.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

tissue_specific_markers: R/tissue_specific_expression.R interesting_gene_reads_wide \
	interesting_isoform_reads_wide
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

interesting_gene_reads: R/subset_reads.R categorized_samples \
	combined_read_counts additional_interesting_genes.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

interesting_isoform_reads: R/subset_reads.R categorized_samples \
	combined_isoform_reads additional_interesting_genes.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

%_reads_wide: %_reads R/wide_reads.R 
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< $@
