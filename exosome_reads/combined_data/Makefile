#!/usr/bin/make -f

### module is how the biocluster loads specific versions; if we're not
### running there, we'll assume the correct version is installed and
### just echo what we're loading
ifdef MODULEPATH
MODULE=module
else
MODULE=echo
endif

SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data

VPATH=R

-include ../srx.mk

SAMPLING=10
READS=1000 5000 10000 50000

GENE_FPKMS:=$(foreach srx,$(SRX_FILES),$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),../$(srx)/$(srx)_split_r$(read)_s$(sample)_genes.fpkm_tracking)))

ISOFORM_FPKMS:=$(foreach srx,$(SRX_FILES),$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),../$(srx)/$(srx)_split_r$(read)_s$(sample)_isoforms.fpkm_tracking)))

ALL_TARGETS=combined_isoform_reads combined_read_counts

all: $(ALL_TARGETS)

combined_read_counts: combine_fpkms.R sample_info_from_name.R $(GENE_FPKMS)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

combined_isoform_reads: combine_isoform_fpkms.R sample_info_from_name.R $(ISOFORM_FPKMS)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

