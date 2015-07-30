#!/usr/bin/make -f

include common_variables.mk
include ../mk/variables.mk

## this file contains all of the SRX and the SRRs in that SRX
include srx.mk

roadmap_epigenomics.tsv:
	wget -O $@ "http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=samples;format=tsv"

human_sequencing_studies.yaml: human_sequencing_studies.pl
	perl $<

human_sequencing_studies.org: yaml_to_org.pl human_sequencing_studies.yaml
	perl $< --include accession --include samples --include description \
		$(wordlist 2,2,$^) > $@

gse_family.mk: human_sequencing_studies.org
	grep ' yes ' $< |awk -F' | ' '{print $$2}'| \
		perl -pe 'BEGIN{print "GSE_FAMILIES="}; s/\n/ /g;'|sed 's/ $$//' \
		> $@

chosen_samples: select_srx.R roadmap_epigenomics.tsv srx_to_sra.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

srx.mk: write_srx.R chosen_samples
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

make_srx_dirs: make_srx_dir.R chosen_samples
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^)

BOWTIE_INDEX_DIR=ref_seqs/

STAR_INDEX_DIR=$(BOWTIE_INDEX_DIR)/star/
star_indexes:
	$(MODULE) load STAR/2.4.2a; \
	mkdir -p $(STAR_INDEX_DIR); \
	STAR --genomeFastaFiles $(BOWTIE_INDEX_DIR)$(FASTA) \
		--runThreadsN $(CORES) \
		--genomeDir $(STAR_INDEX_DIR)

bowtie_indexes: $(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2.1.bt2 $(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2.fa

$(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2.1.bt2: $(BOWTIE_INDEX_DIR)$(FASTA)
	$(MODULE) load bowtie2/2.1.0; \
	bowtie2-build $< $(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2

$(BOWTIE_INDEX_DIR)$(ALIGNMENT_SPECIES)_bt2.fa: $(BOWTIE_INDEX_DIR)$(FASTA)
	rm -f $@
	ln $< $@

$(BOWTIE_INDEX_DIR)$(FASTA): $(BOWTIE_INDEX_DIR)$(FASTA).gz
	gzip -dc $< > $@

$(BOWTIE_INDEX_DIR)$(GTF): $(BOWTIE_INDEX_DIR)$(GTF).gz
	gzip -dc $< > $@

remote_files: $(BOWTIE_INDEX_DIR)$(FASTA).gz $(BOWTIE_INDEX_DIR)$(GTF).gz

fasta: $(BOWTIE_INDEX_DIR)$(FASTA)

gtf: $(BOWTIE_INDEX_DIR)$(GTF)

$(BOWTIE_INDEX_DIR)$(FASTA).gz:
	mkdir -p $(BOWTIE_INDEX_DIR)
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/fasta/$(ALIGNMENT_SPECIES)/dna/$(FASTA).gz" $@

$(BOWTIE_INDEX_DIR)$(GTF).gz:
	mkdir -p $(BOWTIE_INDEX_DIR)
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/gtf/$(ALIGNMENT_SPECIES)/$(GTF).gz" $@


get_srr: $(patsubst %,%-get_srr,$(SRX_FILES))

$(patsubst %,%-get_srr,$(SRX_FILES)): %-get_srr: %
	make -C $* get_srr

.PHONY: get_srr
