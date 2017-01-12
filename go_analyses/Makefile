#!/usr/bin/make -f

VPATH=.:../common_r/
SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl

all: go_tissue_analysis.txt


go_tissue_analysis: go_tissue_analysis.R sample_reads_wide tissue_specific_markers
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

go_tissue_analysis.txt: array_to_text.R go_tissue_analysis
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) go.results $@

sample_reads_wide: tissue_expression.R categorized_samples combined_read_counts
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@


