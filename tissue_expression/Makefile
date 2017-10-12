#!/usr/bin/make -f

include common_variables.mk
include ../mk/variables.mk

## this file contains all of the SRX and the SRRs in that SRX
include srx.mk

## this file contains a variable that describes how to connect to the
## cluster host where I ran the analyses and will retreive them onto
## the local host for analysis
-include cluster_host.mk

roadmap_epigenomics.tsv:
	wget -O $@ "http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=samples;format=tsv"

human_sequencing_studies.yaml: human_sequencing_studies.pl
	perl $<

human_sequencing_studies.org: yaml_to_org.pl human_sequencing_studies.yaml
	perl $< --include accession --include samples --include description \
		$(wordlist 2,2,$^) > $@

GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz:
	wget -O $@ "http://www.gtexportal.org/static/datasets/gtex_analysis_v4/rna_seq_data/$@" || \
		rm $@

GTEx_Analysis_V4_RNA-seq_Flux1.6_transcript_rpkm.txt.gz:
	wget -O $@ "http://www.gtexportal.org/static/datasets/gtex_analysis_v4/rna_seq_data/$@" || \
		rm $@

GTEx_Data_V4_Annotations_SampleAttributesDS.txt:
	wget -O $@ "http://www.gtexportal.org/static/datasets/gtex_analysis_v4/annotations/$@" || \
		rm $@

RME_NOT_IN_SRA=EG.name.txt 57epigenomes.RPKM.pc.gz 57epigenomes.RPKM.nc.gz \
	57epigenomes.RPKM.rb.gz \
	Ensembl_v65.Gencode_v10.ENSG.gene_info

rme_not_in_sra: $(RME_NOT_IN_SRA)

$(RME_NOT_IN_SRA): %: 
	wget -O $@ "http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/$@" || \
		rm $@

# we are currently not selecting the GSE samples this way; we might
# eventually do this later.

# gse_family.mk: human_sequencing_studies.org
# 	grep ' yes ' $< |awk -F' | ' '{print $$2}'| \
# 		perl -pe 'BEGIN{print "GSE_FAMILIES="}; s/\n/ /g;'|sed 's/ $$//' \
# 		> $@

chosen_samples: select_srx.R roadmap_epigenomics.tsv srx_to_sra.txt \
	gse_family_xml/gse_samples
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

srx.mk: write_srx.R chosen_samples
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

make_srx_dirs: make_srx_dir.R chosen_samples
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^)


ISOFORM_FPKMS=$(foreach srx,$(SRX_FILES),$(srx)/$(srx)_isoforms.fpkm_tracking)
GENE_FPKMS=$(foreach srx,$(SRX_FILES),$(srx)/$(srx)_genes.fpkm_tracking)
STAR_LOGS=$(foreach srx,$(SRX_FILES),$(srx)/$(srx)_star/Log.final.out)
combined_read_counts: combine_fpkms.R $(GENE_FPKMS) \
	GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz \
	$(RME_NOT_IN_SRA)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

combined_star_logs: combine_star_logs.R  $(STAR_LOGS)
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

combined_isoform_reads: combine_isoform_fpkms.R $(ISOFORM_FPKMS) \
	GTEx_Analysis_V4_RNA-seq_Flux1.6_transcript_rpkm.txt.gz
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

get_results:
	rsync -avmP \
		--include '**_genes.fpkm_tracking' \
		--include '**_isoforms.fpkm_tracking' \
		--include '**_star/Log.final.out' \
		--include '**_star/ReadsPerGene.out.tab' \
		--include '**_kallisto/abundance.tsv' \
		--include '**_kallisto/run_info.json' \
		--include '**_kallisto.txt' \
		--include '**_kallisto_info.txt' \
		--include '**trinity_diamond.txt' \
		--include '**_trinity_align_rsem_isoforms.txt' \
        --include '**_fastqc.html' \
        --include '**_fastqc.zip' \
        --include '**/' \
	    --exclude '**' \
		$(CLUSTER_HOST)/ .;

submit_trimmed_fastqc: $(patsubst %,%-submit_trimmed_fastqc,$(SRX_FILES))

submit_alignment: $(patsubst %,%-submit_alignment,$(SRX_FILES))

submit_call: $(patsubst %,%-submit_call,$(SRX_FILES))

submit_kallisto_call: $(patsubst %,%-submit_kallisto_call,$(SRX_FILES))

$(patsubst %,%-submit_trimmed_fastqc,$(SRX_FILES)): %-submit_trimmed_fastqc:
	+make -C $* submit_trimmed_fastqc

$(patsubst %,%-submit_alignment,$(SRX_FILES)): %-submit_alignment:
	+make -C $* submit_alignment

$(patsubst %,%-submit_call,$(SRX_FILES)): %-submit_call:
	+make -C $* submit_call

$(patsubst %,%-submit_kallisto_call,$(SRX_FILES)): %-submit_kallisto_call:
	+make -C $* submit_kallisto_call

get_srr: $(patsubst %,%-get_srr,$(SRX_FILES))

$(patsubst %,%-get_srr,$(SRX_FILES)): %-get_srr: %
	make -C $* get_srr

.PHONY: get_srr
