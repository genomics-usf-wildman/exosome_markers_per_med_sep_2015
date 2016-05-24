#!/usr/bin/perl
# read_biaser.pl biases reads towards motifs
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2016 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

read_biaser.pl - biases reads towards motifs

=head1 SYNOPSIS

read_biaser.pl [options]

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

read_biaser.pl

=cut


use IO::Uncompress::AnyUncompress;
use vars qw($DEBUG);

# In T cells, GGAG is required to target miRNA to T cell
# motifs\cite{Villarroya-Beltri.ea2013:SumoylatedhnRNPA2B1controlssorting}.

# ACCAGCCU, CAGUGAGC, and UAAUCCCA motifs account for about 25% of
# eRNAs. \cite{Batagov.ea2013:Exosomessecretedbyhuman}


my %options = (debug           => 0,
               help            => 0,
               man             => 0,
               seed            => rand,
               reads           => undef,
               skip_reads      => 0,
               count           => 0,
               min_likelihood => 0.5,
               motif_output_abundance => 0.25,
               motif_input_abundance  => 0.01,
               motif => [qw(ACCAGCCU CAGUGAGC UAAUCCCA)],
              );

GetOptions(\%options,
           'reads=s',
           'seed|random-seed=s',
           'count!',
           'skip_reads|skip-reads=i',
           'motif_input_abundance|motif-input-abundance=s',
           'motif_output_abundance|motif-output-abundance=s',
           'motif=s@',
           'output=s',
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;
if (defined $options{reads} and $options{reads} < 1) {
    push @USAGE_ERRORS,"You should give at least one read";
}

srand($options{seed});

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;


my @files = @ARGV;
if (not @ARGV) {
    @files = undef;
}

# build motif markers regex
my $motif_regex = '('.join('|',map {s/[UT]/\[TU\]/; $_} @{$options{motif}}).')';
print STDERR "motif regex: $motif_regex\n";
$motif_regex = qr/$motif_regex/;

my $output = \*STDOUT;
if (defined $options{output}) {
    open($output,'>',$options{output}) or
        die "Unable to open $options{output} for writing: $!";
}

# in order to have output abundance, motif containing reads need to
# have more (or less) probability of being output than normal reads.

my $motif_keep = $options{min_likelihood};
my $non_motif_keep = $options{min_likelihood};
if ($options{motif_output_abundance} > $options{motif_input_abundance}) {
    $non_motif_keep = $motif_keep / ($options{motif_output_abundance} /
                                     $options{motif_input_abundance});
} else {
    $motif_keep = $non_motif_keep / ((1-$options{motif_output_abundance}) /
                                     (1-$options{motif_input_abundance}));
}

my $count = 0;
my $total_reads_kept = 0;
for my $fn (@files) {
    my $fh;
    if (not defined $fn) {
        $fh = \*STDIN;
    } else {
        $fh = IO::Uncompress::AnyUncompress->new($fn,MultiStream => 1) or
            die "Unable to open $fn for reading: $!";
    }
    my $read_pos = 9999;
    my $current_read = '';
    my $keep_read = 0;
    while (<$fh>) {
        $read_pos++;
        if ($read_pos > 3) {
            $read_pos = 0;
        }
        if ($read_pos == 0) {
            $current_read = $_;
            $keep_read = 0;
        } elsif ($read_pos == 1) {
            if ($options{count}) {
                if ($_ =~ $motif_regex) {
                    $count++;
                }
                $total_reads_kept++;
                last if defined $options{reads} and $total_reads_kept > $options{reads};
            } else {
                my $r = rand;
                # first, we want to ignore a certain percentage of all reads
                if ($r <= $options{min_likelihood}) {
                    # then, if we're going to keep it at all, check to see
                    # if it matches a motif.
                    if ($_ =~ $motif_regex) {
                        print STDERR "matched motif; ";
                        # it matches;
                        if ($r <= $motif_keep) {
                            print STDERR "kept ";
                            $keep_read = 1;
                        }
                        print STDERR "($r)\n";
                    } else {
                        print STDERR "didn't match motif; ";
                        if ($r <= $non_motif_keep) {
                            print STDERR "kept ";
                            $keep_read = 1;
                        }
                        print STDERR "($r)\n";
                    }
                }
            }
        }
        if ($keep_read) {
            $current_read .= $_;
            if ($read_pos == 3) {
                print {$output} $current_read;
                $total_reads_kept++;
                if (defined $options{reads} and $total_reads_kept >= $options{reads}) {
                    last;
                }
            }
        }
    }
}

if ($options{count}) {
    print {$output} "motif_input_abundance = ".$count/$total_reads_kept;
}

close($output);




__END__
