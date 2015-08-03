#!/usr/bin/make -f

BOWTIE_OPTIONS?=
BOWTIE_INDEX_DIR?=../ref_seqs/
STAR_INDEX_DIR?=$(BOWTIE_INDEX_DIR)/star/
LOCAL_BOWTIE_OPTIONS?= -N 1 -L 15 -k 10 --local
ALIGNMENT_SPECIES?=$(SPECIES)

# we need to use a comma in a rule below, so this handles that
# escaping.
comma=,

call: $(SRX)_genes.fpkm_tracking

$(SRX)_genes.fpkm_tracking: $(SRX)_star.bam $(BOWTIE_INDEX_DIR)$(GTF)
	$(MODULE) load cufflinks/2.2.1; \
	cufflinks -p $(CORES) -G $(wordlist 2,2,$^) $<
	for file in genes.fpkm_tracking isoforms.fpkm_tracking skipped.gtf transcripts.gtf; do \
		mv $${file} $(SRX)_$${file}; \
	done;

alignment: $(SRX)_star.bam

$(SRX)_star.bam:
	$(MODULE) load STAR/2.4.2a; \
	mkdir -p $(SRX)_star; \
	STAR --outFileNamePrefix $(SRX)_star/ \
		--outSAMtype BAM SortedByCoordinate \
		--runThreadN $(CORES) \
        --outSAMstrandField intronMotif \
		--genomeDir $(STAR_INDEX_DIR) \
		--readFilesCommand "gzip -c" \
		--readFilesIn $(TOPHAT_FASTQ_ARGUMENT)

$(SRX)_tophat.bam: \
	$(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2.1.bt2 \
	$(FASTQ_FILES) $(BOWTIE_INDEX_DIR)$(GTF) $(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2.fa
	$(MODULE) load tophat2/2.0.10; \
	tophat -G $(BOWTIE_INDEX_DIR)$(GTF) -p $(CORES) \
		--transcriptome-index $(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_tophat_indexes \
		-o $(SRX)_tophat \
		$(patsubst %.1.bt2,%,$(wordlist 1,1,$^)) \
		$(TOPHAT_FASTQ_ARGUMENT);
	ln $(SRX)_tophat/accepted_hits.bam $@ -s

