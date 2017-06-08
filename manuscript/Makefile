#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

MANUSCRIPT_TEX_DEPS=../tissue_specific_expression/tissue_specific_markers	\
../tissue_specific_expression/interesting_gene_reads_wide					\
../tissue_specific_expression/combined_star_logs							\
../tissue_specific_expression/svm_test_results								\
../tissue_specific_expression/svm_test_results_reduced						\
../tissue_specific_expression/housekeeping_genes_superset


MANUSCRIPT_PDF_DEPS=acronyms.tex

all: manuscript.pdf

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(MANUSCRIPT_PDF_DEPS)
	latexmk -pdf -f -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

manuscript.tex: manuscript.Rnw $(MANUSCRIPT_TEX_DEPS)

