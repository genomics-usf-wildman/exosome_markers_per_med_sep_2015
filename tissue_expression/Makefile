#!/usr/bin/make -f

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
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

srx.mk: write_srx.R chosen_samples
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

make_srx_dirs: make_srx_dir.R chosen_samples
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^)
