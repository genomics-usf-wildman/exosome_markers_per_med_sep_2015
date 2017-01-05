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
else
SRR_FASTQ_FILES=$(foreach srr,$(SRRS),$(srr)_1.fastq.gz $(srr)_2.fastq.gz)
endif

make_srr_fastq: $(SRR_FASTQ_FILES)

ifeq ($(NREADS),1)
$(SRR_FASTQ_FILES): %.fastq.gz:
	$(MODULE) load sratoolkit/2.3.5-2; \
	fastq-dump --split-3 -B --gzip $*;
else
%_1.fastq.gz %_2.fastq.gz:
	$(MODULE) load sratoolkit/2.3.5-2; \
	fastq-dump --split-3 -B --gzip $*;
endif

SAMPLING=10
READS=1000 5000 10000 50000
FASTQ_FILES=$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),$(SRX)_r$(read)_s$(sample).fastq.gz))

STAR_OPTIONS=--genomeLoad LoadAndKeep --sjdbGTFfile $(GTF) --quantMode GeneCounts

ifeq ($(NREADS),1)
READ_BIASER_OPTS:=
else
READ_BIASER_OPTS:=--paired
endif


make_fastq: ../read_biaser.pl $(SRR_FASTQ_FILES)
	$(MODULE) load perl/5.20.1; \
	$< $(foreach read,$(READS),--read $(read)) --samplings $(SAMPLING) \
		--output-prefix $(SRX) $(READ_BIASER_OPTS) $(SRR_FASTQ_FILES)
	gzip *_r*_s*.fastq;
	touch $@

include ../../rnaseq_workflow/common_makefile
