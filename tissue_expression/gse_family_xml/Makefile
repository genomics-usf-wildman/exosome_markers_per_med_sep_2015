#!/usr/bin/make -f

SHELL=/bin/bash

VPATH=../../common_r

## GSE_FAMILIES
include ../gse_family.mk
R=R
ROPTS=-q --no-save --no-restore-data

GSE_FAMILIES_TXT=$(patsubst %,%_family.txt,$(GSE_FAMILIES))

GSE_FAMILIES_XML=$(patsubst %,%_family.xml,$(GSE_FAMILIES))
GSE_FAMILIES_XMLTGZ=$(patsubst %,%_family.xml.tgz,$(GSE_FAMILIES))

gse_samples.txt: array_to_text.R gse_samples 
	$(R) $(ROPTS) -f $< --args gse_samples gse.samples $@

gse_samples: combine_gse_families.R $(GSE_FAMILIES_TXT)
	$(R) $(ROPTS) -f $< --args $(GSE_FAMILIES_TXT) $@

gse_families_txt: $(GSE_FAMILIES_TXT)

$(GSE_FAMILIES_TXT): %_family.txt: %_family.xml ../parse_family_xml.pl
	../parse_family_xml.pl $< > $@ 

gse_families_xml: $(GSE_FAMILIES_XML)

$(GSE_FAMILIES_XML): %_family.xml: %_family.xml.tgz
	tar -zmxf $< $@

%_family.xml.tgz:
	wget -c "ftp://ftp.ncbi.nlm.nih.gov/geo/series/$$(echo -n $*|perl -pe 's/\d{0,3}$$/nnn/')/$*/miniml/$*_family.xml.tgz"
