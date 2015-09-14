VPATH=R:perl
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL=perl

ifdef MODULEPATH
MODULE=module
else
MODULE=echo
endif
