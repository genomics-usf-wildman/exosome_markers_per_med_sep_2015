SPECIES=homo_sapiens
ALIGNMENT_SPECIES=homo_sapiens
GTF=Homo_sapiens.GRCh38.80.gtf
FASTA=Homo_sapiens.GRCh38.dna.toplevel.fa
SHELL=/bin/bash

### module is how the biocluster loads specific versions; if we're not
### running there, we'll assume the correct version is installed and
### just echo what we're loading
ifdef MODULEPATH
MODULE=module
else
MODULE=echo
endif

# this is the env variable to tell us how many processors on this node
# we get
ifdef PBS_NUM_PPN
CORES=$(PBS_NUM_PPN)
else
CORES=8
endif
