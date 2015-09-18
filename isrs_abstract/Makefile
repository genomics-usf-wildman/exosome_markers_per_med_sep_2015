#!/usr/bin/make -f

all: isrs_abstract.doc

EMACSARGS=--batch  -l ~/.emacs 
EMACS=emacs

ORG_FILES?=$(wildcard *.org)

PDF_FROM_ORG_FILES=$(patsubst %.org,%.pdf,$(ORG_FILES))

DOC_FROM_ORG_FILES=$(patsubst %.org,%.doc,$(ORG_FILES))

$(DOC_FROM_ORG_FILES): %.doc: %.org
	$(EMACS) $(EMACSARGS) --visit $< --eval '(setq org-confirm-babel-evaluate nil)' \
		--eval '(setq org-odt-preferred-output-format "doc")' \
		--eval '(setq ess-ask-for-ess-directory nil)' \
		--funcall org-odt-export-to-odt
