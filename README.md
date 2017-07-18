# Tissue Specific Markers in Extracellular Vesicles

This repository contains code which identifies tissue specific markers
based on the GTEx and roadmap epigenomics project databases, as well
as additional data from SRA and GEO.

The primary question that is being asked (and hopefully answered) is:

1. How can we identify the source tissue of exosomes found in a
   peripheral fluid?
    2. What are genes whose expression is particular to a tissue?
    3. Knowing those genes, can we identify the source tissue of *in
       silico* generated exosome reads?

# Contents of the Repo

1. `common_r` contains various R utilities for sub analyses
2. `exosome_reads` contains code to download placenta and various
   other tissue SRA experiments, perform an alignment, bias and subset
   the aligned reads to exosomes, and then call the resultant
   subset+biased bam files using cufflinks, and then pile all of the
   data up. See its `Makefile` for details.
3. `exosome_sra` contains code to examine reads from various
   extracellular vesicles from different tissue types. Currently not
   used for anything.
4. `go_analyses` contains the code to perform GO analyses on the basis
   of the identified tissue specific markers. See its `Makefile` for
   details.
5. `isrs_abstract` is an abstract for the Illinois Symposium on
   Reproductive Science which was accepted
6. `manuscript` contains the code to produce the paper based upon this
   work. See its `Makefile` for details.
7. `mk` contains various makefile settings
8. `poster` contains the code for a poster which was presented at the
   Mayo Individualizing Medicine conference.
9. `rnaseq_workflow` is a git submodule which is the standard RNAseq
   workflow
10. `slides` is a selection of slides for the ISRS talk
11. `tissue_expression` contains all of the code to build an R matrix
    of tissue expression across a large number of samples. See its
    `Makefile` for details.
12. `tissue_specific_expression` contains the code to turn the matrix
    of tissue expression into a list of genes which are possibly
    tissue specific and genes which are not by calculating τ. It also
    contains code to run SVM and other routines using `caret` to do
    machine learning to identify source tissues.

