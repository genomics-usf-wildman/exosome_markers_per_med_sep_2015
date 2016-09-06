#!/usr/bin/make -f

include ../mk/variables.mk

include ../tissue_expression/common_variables.mk

include ../tissue_expression/srx.mk

-include cluster_host.mk

make_srx_dirs: make_srx_dir.R ../tissue_expression/chosen_samples
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^)

make_fastq: $(patsubst %,%-make_fastq,$(SRX_FILES))

$(patsubst %,%-make_fastq,$(SRX_FILES)): %-make_fastq: %
	make -C $* make_fastq
