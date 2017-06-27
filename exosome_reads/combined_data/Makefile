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

NORMALIZE_FILES=combined_read_counts combined_isoform_reads

NORMALIZED_DATA_FILES=$(foreach file,$(NORMALIZE_FILES),$(file)_normalized)

ALL_TARGETS=hugo_nomenclature $(NORMALIZED_DATA_FILES)	\
$(NORMALIZE_FILES)

all: $(ALL_TARGETS)

combined_read_counts: combine_fpkms.R sample_info_from_name.R $(GENE_FPKMS)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

combined_isoform_reads: combine_isoform_fpkms.R sample_info_from_name.R $(ISOFORM_FPKMS)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

combined_fastqc: combine_fastqc.R sample_info_from_name.R $(FASTQC_FILES)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

hugo_nomenclature.txt:
	wget -O $@ "http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_ensembl_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"

hugo_nomenclature: hugo_nomenclature.R hugo_nomenclature.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

$(NORMALIZED_DATA_FILES): %_normalized: % quantile_normalize.R
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $< $@
