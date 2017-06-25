#!/usr/bin/make -f

SPECIES=homo_sapiens
ALIGNMENT_SPECIES=homo_sapiens
ENSEMBL_RELEASE=80
SHELL=/bin/bash

STRIP_PATCHES:=1
STRIP_PATCHES_SCRIPT:=../../rnaseq_workflow/strip_patches.pl
STRIP_PATCHES_OPTIONS:=--valid-chr='^(?:chr)?[12]?\d|MT|[XY]$$'
REFERENCE_DIR:=../../tissue_expression/ref_seqs/
NOTRIM=1
include srx_info.mk

ifeq ($(NREADS),1)
SRR_FASTQ_FILES=$(foreach srr,$(SRRS),$(srr).fastq.gz)
SRX_FASTQ_FILES=$(SRX).fastq.gz
else
SRR_FASTQ_FILES=$(foreach srr,$(SRRS),$(srr)_1.fastq.gz $(srr)_2.fastq.gz)
SRX_FASTQ_FILES=$(SRX)_1.fastq.gz $(SRX)_2.fastq.gz
endif

make_srr_fastq: $(SRR_FASTQ_FILES)

make_srx_fastq: $(SRX_FASTQ_FILES)

ifeq ($(NREADS),1)
$(SRR_FASTQ_FILES): %.fastq.gz:
	$(MODULE) load sratoolkit/2.3.5-2; \
	fastq-dump --split-3 -B --gzip $*;
else
%_1.fastq.gz %_2.fastq.gz:
	$(MODULE) load sratoolkit/2.3.5-2; \
	fastq-dump --split-3 -B --gzip $*;
endif

ifeq ($(NREADS),1)
$(SRX).fastq.gz: $(SRR_FASTQ_FILES)
	pigz -dc $^ |pigz -c > $@
else
$(SRX)_1.fastq.gz: $(foreach srr,$(SRRS),$(srr)_1.fastq.gz)
	pigz -dc $^ |pigz -c > $@

$(SRX)_2.fastq.gz: $(foreach srr,$(SRRS),$(srr)_2.fastq.gz)
	pigz -dc $^ |pigz -c > $@
endif



SAMPLING=10
READS=1000 5000 10000 50000
FASTQ_FILES:=$(SRX_FASTQ_FILES)

# we'd like to use shared memory, but that's not supported with a GTFfile
STAR_OPTIONS=--sjdbGTFfile $(GTF) --quantMode GeneCounts

# make_fastq: ../read_biaser.pl $(SRR_FASTQ_FILES)
# 	$(MODULE) load perl/5.20.1; \
# 	$< $(foreach read,$(READS),--read $(read)) --samplings $(SAMPLING) \
# 		--output-prefix $(SRX) $(READ_BIASER_OPTS) $(SRR_FASTQ_FILES)
# 	gzip *_r*_s*.fastq;
# 	touch $@

FPKM_GENES_ANALYSIS_FILES:=

SPLIT_FPKM_GENES_ANALYSIS_FILES:=$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),$(SRX)_split_r$(read)_s$(sample)_genes.fpkm_tracking))

split_call: $(SPLIT_FPKM_GENES_ANALYSIS_FILES)

%_genes.fpkm_tracking %_isoforms.fpkm_tracking %_skipped.gtf %_transcripts.gtf: %.bam \
	$(GTF)
	mkdir -p $(*)_cufflinks;
	$(MODULE) load cufflinks/2.2.1; \
	cufflinks -o $(*)_cufflinks $(CUFFLINKS_OPTIONS) -p $(CORES) -G $(wordlist 2,2,$^) $<
	for file in genes.fpkm_tracking isoforms.fpkm_tracking skipped.gtf transcripts.gtf; do \
		mv $(*)_cufflinks/$${file} $(*)_$${file}; \
	done;
	rm $(*)_cufflinks -rf;

include ../../rnaseq_workflow/common_makefile

SPLIT_STAR_ALIGNMENT_FILES:=$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),$(SRX)_split_r$(read)_s$(sample).bam))

SPLIT_STAR_ALIGNMENT_FILES_PATTERN:=$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),%_split_r$(read)_s$(sample).bam))

split_bams: $(SPLIT_STAR_ALIGNMENT_FILES)

$(SPLIT_STAR_ALIGNMENT_FILES_PATTERN): %_star.bam ../read_biaser_bam.pl 
	$(MODULE) load samtools; \
	$(MODULE) load perl; \
	../read_biaser_bam.pl $(foreach read,$(READS),--read $(read)) --samplings $(SAMPLING) \
		--output-prefix $(SRX)_split $(STAR_ALIGNMENT_FILES)

