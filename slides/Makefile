#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: slides.pdf

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(wildcard *.tex)
	latexmk -f -pdf -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

slides.tex: slides.Rnw
