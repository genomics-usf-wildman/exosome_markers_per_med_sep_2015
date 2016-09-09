#!/usr/bin/make -f

SPECIES=homo_sapiens
ALIGNMENT_SPECIES=homo_sapiens
ENSEMBL_RELEASE=80
SHELL=/bin/bash

STRIP_PATCHES:=1
STRIP_PATCHES_SCRIPT:=../../rnaseq_workflow/strip_patches.pl
STRIP_PATCHES_OPTIONS:=--valid-chr='^(?:chr)?[12]?\d|MT|[XY]$$'
REFERENCE_DIR:=../../tissue_expression/ref_seqs/
include srx_info.mk


SRR_FASTQ_FILES=$(patsubst %, %.fastq.gz,$(SRRS))

make_srr_fastq: $(SRR_FASTQ_FILES)

$(SRR_FASTQ_FILES): %.fastq.gz:
	$(MODULE) load sratoolkit/2.3.5-2; \
	fastq-dump -B --gzip $*;


SAMPLING=10
READS=1000 5000 10000 50000
FASTQ_FILES=$(foreach sample,$(shell seq 1 $(SAMPLING)),$(foreach read,$(READS),$(SRX)_r$(read)_s$(sample).fastq.gz))

make_fastq: ../read_biaser.pl $(SRR_FASTQ_FILES)
	$< $(foreach read,$(READS),--read $(read)) --samplings $(SAMPLING) \
		--output-prefix $(SRX) $(SRR_FASTQ_FILES)
	gzip *_r*_s*.fastq;
	touch $@

include ../../rnaseq_workflow/common_makefile
