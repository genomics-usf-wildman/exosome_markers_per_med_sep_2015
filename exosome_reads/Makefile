#!/usr/bin/make -f

include ../mk/variables.mk

include ../tissue_expression/common_variables.mk

include srx.mk

-include cluster_host.mk

all:

make_srx_dirs: make_srx_dir.R ../tissue_expression/chosen_samples		\
../tissue_specific_expression/categorized_samples chosen_tissues.yaml
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^)

make_fastq: $(patsubst %,%-make_fastq,$(SRX_FILES))

make_srr_fastq: $(patsubst %,%-make_srr_fastq,$(SRX_FILES))

$(patsubst %,%-make_fastq,$(SRX_FILES)): %-make_fastq: %
	+make -C $* make_fastq

$(patsubst %,%-make_srr_fastq,$(SRX_FILES)): %-make_srr_fastq: %
	+make -C $* make_srr_fastq

srx.mk: write_srx.R ../tissue_expression/chosen_samples		\
../tissue_specific_expression/categorized_samples chosen_tissues.yaml
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

