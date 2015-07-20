#!/usr/bin/make

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
