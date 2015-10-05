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


brain_specific_markers: R/brain_specific.R categorized_samples combined_read_counts \
	brain_tissues.yaml \
	housekeeping_genes
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

tissue_specific_markers: R/tissue_specific_expression.R interesting_gene_reads_wide \
	interesting_isoform_reads_wide
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

housekeeping_genes: R/housekeeping_genes.R categorized_samples combined_read_counts
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

tissue_specific_markers.txt: tissue_specific_markers ../common_r/array_to_text.R
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< tissue.specific.genes $@

brain_specific_markers.txt: brain_specific_markers ../common_r/array_to_text.R
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< brain.tissue.specific.genes $@

housekeeping_genes.txt: housekeeping_genes ../common_r/array_to_text.R
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< housekeeping.genes $@

interesting_gene_reads: R/subset_reads.R categorized_samples \
	combined_read_counts additional_interesting_genes.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

interesting_isoform_reads: R/subset_reads.R categorized_samples \
	combined_isoform_reads additional_interesting_genes.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

housekeeping_genes_eisenberg_2003.html:
	wget -O $@ "https://web.archive.org/web/20070805203436/http://compugen.co.il/supp_info/Housekeeping_genes.html"

.ONESHELL:
housekeeping_genes_eisenberg_2003.txt: housekeeping_genes_eisenberg_2003.html
	perl -0777 -ne 'while (m/\*?(NM_\d+)<[^>]+>(?:\&nbsp\;)*\s*<[^>]+>\s*([^\<]*)</msg) { 
		my $$nm = $$1;
		my $$gene = $$2;
	 	$$nm =~ s/[\r\n]+/ /g;
	 	$$gene =~ s/[\r\n]+/ /g;
		my ($$name,$$long) = $$gene =~ m/(.+)\s+\(([^\)]+)\),/; 
	    print qq($$nm\n); }' $< > $@

housekeeping_genes_eisenberg_2013.txt:
	wget -O $@ "http://www.tau.ac.il/~elieis/HKG/HK_genes.txt"

housekeeping_genes_eisenberg_2013_annotated.txt: R/annotate_nm.R housekeeping_genes_eisenberg_2013.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

housekeeping_genes_eisenberg_2003_annotated.txt: R/annotate_nm.R housekeeping_genes_eisenberg_2003.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

housekeeping_genes_superset: R/housekeeping_genes_superset.R housekeeping_genes housekeeping_genes_eisenberg_2013_annotated.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

housekeeping_genes_superset.txt: housekeeping_genes_superset ../common_r/array_to_text.R
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< housekeeping.genes.superset $@



%_reads_wide: %_reads R/wide_reads.R 
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< $@
