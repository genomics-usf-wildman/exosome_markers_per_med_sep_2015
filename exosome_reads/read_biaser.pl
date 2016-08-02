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

use List::Util qw(sum);
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
               min_likelihood => 0.5,
               motif_output_abundance => 0.25,
               # motif_input_abundance  => 0.01,
               motif => [qw(ACCAGCCU CAGUGAGC UAAUCCCA)],
              );

GetOptions(\%options,
           'reads=s',
           'seed|random-seed=s',
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

my $count = 0;
my $total_reads_kept = 0;
my @file_counts;
my @total_counts;
use Data::Printer;
my @file_reads;
for my $fn (@files) {
    my $fh;
    if (not defined $fn) {
        $fh = \*STDIN;
    } else {
        p $fn;
        $fh = IO::Uncompress::AnyUncompress->new($fn,MultiStream => 1) or
            die "Unable to open $fn for reading: $!";
    }
    my $read_pos = 9999;
    my $current_read;
    my $read_id;
    my $match_motif;
    my $quality_ok = 0;
    push @file_reads,[];
    push @file_counts,[0,0];
    while (<$fh>) {
        $read_pos++;
        if ($read_pos > 3) {
            $read_pos = 0;
            if ($quality_ok and defined $read_id) {
                if ($match_motif) {
                    # store matching read ids and count them per file
                    push @{$file_reads[-1][0]},$read_id;
                    $file_counts[-1][0]++;
                    $total_counts[0]++;
                } else {
                    push @{$file_reads[-1][1]},$read_id;
                    $file_counts[-1][1]++;
                    $total_counts[1]++;
                }
            }
        }
        if ($read_pos == 0) {
            $current_read = $_;
            ($read_id) = $_ =~ /^(\S+)/;
        } elsif ($read_pos == 1) {
            $match_motif = $_ =~ $motif_regex;
        } elsif ($read_pos == 3) {
            # if the average quality is more than 15, keep the read.
            # [That's pretty generous, actually.]
            $quality_ok = sum(map {ord($_) - 33} split '',$_)/length($_) > 15;
        }

    }
    close($fh);
}

my @file_keep_reads;
## figure out the total reads and the proportion of matched and non-matched
print STDERR "There were ".sum(@total_counts)." reads in ".scalar @file_counts." files\n";
print STDERR $total_counts[0]." contained the motif, ".$total_counts[1]." did not\n";
my $kept_reads = 0;
while ($kept_reads < $options{reads}) {
    my $m = rand > $options{motif_output_abundance} ? 1 : 0;
    # choose a motif matching read ($m==0) or non-matching ($m==1)
    my $r = rand $total_counts[$m];
    my $f = 0;
    while ($r > $file_counts[$f][$m]) {
        $r -= $file_counts[$f][$m];
    }
    if ($file_keep_reads[$f]{$file_reads[$f][$m][$r]}) {
        next;
    } else {
        $file_keep_reads[$f]{$file_reads[$f][$m][$r]} = 1;
        $kept_reads++;
    }
}

my $f = 0;
for my $fn (@files) {
    my $fh;
    if (not defined $fn) {
        $fh = \*STDIN;
    } else {
        $fh = IO::Uncompress::AnyUncompress->new($fn,MultiStream => 1) or
            die "Unable to open $fn for reading: $!";
    }
    my $read_pos = 999;
    my $keep_read = 0;
    while (<$fh>) {
        $read_pos++;
        if ($read_pos > 3) {
            $read_pos = 0;
        }
        if ($read_pos == 0) {
            my ($read_id) = $_ =~ /^(\S+)/;
            $keep_read = exists $file_keep_reads[$f]{$read_id};
        }
        print {$output} $_ if $keep_read;
    }
    close($fh);
    $f++;
}





__END__
