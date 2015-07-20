#!/usr/bin/perl
# human_sequencing_studies.pl retreives and parses information on human sequencing studies
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2014 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

human_sequencing_studies.pl -  retreives and parses information on human sequencing studies

=head1 SYNOPSIS

human_sequencing_studies.pl [options]

 Options:
   --debug, -d debugging level (Default 0)
   --help, -h display this help
   --man, -m display manual

=head1 OPTIONS

=over

=item B<--debug, -d>

Debug verbosity. (Default 0)

=item B<--help, -h>

Display brief usage information.

=item B<--man, -m>

Display this manual.

=back

=head1 EXAMPLES

human_sequencing_studies.pl

=cut


use vars qw($DEBUG);

use WWW::Mechanize;
use IO::File;
use YAML;

my %options = (debug           => 0,
               help            => 0,
               man             => 0,
              );

GetOptions(\%options,
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;
# if (1) {
#     push @USAGE_ERRORS,"You must pass something";
# }

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;

sub get_experiments {
    # temporarily use this for debugging
#     my $fh = IO::File->new('human_sequencing_studies.txt','<:encoding(UTF-8)') or
#         die "Unable to open human_sequencing_studies.txt for reading: $!";
#     my $ret = '';
#     while (<$fh>) {
#         $ret .= $_;
#     }
#     close $fh;
#     return $ret;
    my $m = WWW::Mechanize->new();

    $m->get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE[ETYP]+AND+human[Organism]+AND+expression+profiling+by+high+throughput+sequencing[DataSet+Type]&usehistory=y');
    my ($queryid) = $m->content() =~ /<QueryKey>(\d+)\<\/QueryKey>/i;
    my ($webenv) = $m->content() =~ /<WebEnv>([^<]+)\<\/WebEnv>/i;
    $m->get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&query_key='.
            $queryid.
            '&WebEnv='.$webenv);

    return($m->content());
}

my %experiments;
my %current_experiment;

for my $line (split /\n/,get_experiments()) {
    local $_ = $line;
    if (/^$/) {
        $experiments{$current_experiment{accession}} = {%current_experiment} if
            defined $current_experiment{accession};
        %current_experiment = ();
        next;
    }
    if (/^(\d+)\.\s+(.+)$/) {
        $current_experiment{number} = $1;
        $current_experiment{description} = $2;
    } elsif (/\(Submitter\s+supplied\)\s+(.+)$/) {
        $current_experiment{description_extra} = $1;
    } elsif (/^Organism:\s+(.+)$/) {
        $current_experiment{organism} = $1;
    } elsif (/^Series\s+Accession:\s+([^\s]+)\s+ID:\s+(\d+)\s*$/) {
        $current_experiment{accession} = $1;
        $current_experiment{id} = $2;
    } elsif (/^Platforms?:\s+(.+?)\s+(\d+)\s+Samples?\s*$/) {
        $current_experiment{samples} = $2;
        $current_experiment{platform} = $1;
    } elsif (/(\d+)\s+Samples?\s*$/) {
        $current_experiment{samples} = $1;
    }
}
$experiments{$current_experiment{accession}} = {%current_experiment} if
    defined $current_experiment{accession};

YAML::DumpFile('human_sequencing_studies.yaml',
               [map {$experiments{$_}}
                sort {($experiments{$b}{samples}//0) <=> ($experiments{$a}{samples}//0)}
                keys %experiments]);

__END__
