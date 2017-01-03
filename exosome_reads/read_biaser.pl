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

read_biaser.pl - biases reads towards motifs by sampling reads without replacement

=head1 SYNOPSIS

read_biaser.pl [options]

 Options:
   --reads number of reads to sample
   --output-prefix name of output file prefix
   --seed,--random-seed  Random seed
   --samplings Number of times to sample reads
   --motif-output-abundance Abundance of motif in output (0.25)
   --debug, -d debugging level (Default 0)
   --help, -h display this help
   --man, -m display manual

=head1 OPTIONS

=over

=item B<--reads>

Number of reads to sample. If repeated, each read option will be
output C<samplings> times. Required.

=item B<--samplings>

Number of samplings to perform; defaults to 1. [Can be used to sample
multiple sets of reads from a single file without re-reading the
file.]

=item B<--output-prefix>

Output prefix of file. Defaults to
biased_reads_rC<reads>_sC<sample>.fastq where C<reads> is the number
of reads, and C<sample> is the nth sampling.

=item B<--seed,--random-seed>

Random seed to use; defaults to rand()

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
               samplings       => 1,
               reads           => [],
               skip_reads      => 0,
               min_likelihood => 0.5,
               'output_prefix' => 'biased_reads',
               motif_output_abundance => 0.25,
               # motif_input_abundance  => 0.01,
               motif => [qw(ACCAGCCU CAGUGAGC UAAUCCCA)],
              );

GetOptions(\%options,
           'reads=i@',
           'seed|random-seed=s',
           'motif_output_abundance|motif-output-abundance=s',
           'motif=s@',
           'samplings=i',
           'output_prefix|output-prefix=s',
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;
if (not @{$options{reads}} or
    grep {$_ <= 1} @{$options{reads}}
   ) {
    push @USAGE_ERRORS,"You must pass --reads and all --reads must be >= 1";
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

my @output_files;
for my $i (0..$#{$options{reads}}) {
    for my $j (0..($options{samplings}-1)) {
        my $fn = $options{output_prefix}.
            '_r'.$options{reads}[$i].
            '_s'.($j+1).'.fastq';
        $output_files[$i][$j] =
            IO::File->new($fn,'w') or
                die "Unable to open $fn for writing: $!";
    }
}


## read in the file to figure out which reads match the motifs, and
## which do not
my $count = 0;
my $total_reads_kept = 0;
my @file_counts;
my @total_counts;
my @file_reads;
for my $fn (@files) {
    my $fh;
    if (not defined $fn) {
        $fh = \*STDIN;
    } else {
        open($fh,'-|','gzip','-dc',$fn) or
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
for my $i (0..$#{$options{reads}}) {
    my $read_count = $options{reads}[$i];
    for my $j (0..($options{samplings}-1)) {
        my $kept_reads = 0;
        while ($kept_reads < $read_count) {
            my $m = rand > $options{motif_output_abundance} ? 1 : 0;
            # choose a motif matching read ($m==0) or non-matching ($m==1)
            # at random from the total set of reads
            my $r = rand $total_counts[$m];
            my $f = 0;
            while ($r > $file_counts[$f][$m]) {
                $r -= $file_counts[$f][$m];
            }
            if ($file_keep_reads[$f]{$file_reads[$f][$m][$r]}{$i.'_'.$j}) {
                next;
            } else {
                $file_keep_reads[$f]{$file_reads[$f][$m][$r]}{$i.'_'.$j} = 1;
                $kept_reads++;
            }
        }
    }
}

my $f = 0;
for my $fn (@files) {
    my $fh;
    if (not defined $fn) {
        $fh = \*STDIN;
    } else {
        open($fh,'-|','gzip','-dc',$fn) or
            die "Unable to open $fn for reading: $!";
    }
    my $read_pos = 999;
    my $keep_read_id = 0;
    while (<$fh>) {
        $read_pos++;
        if ($read_pos > 3) {
            $read_pos = 0;
        }
        if ($read_pos == 0) {
            $keep_read_id = 0;
            my ($read_id) = $_ =~ /^(\S+)/;
            if (exists $file_keep_reads[$f]{$read_id}) {
                $keep_read_id = $read_id;
            }
        }
        if ($keep_read_id) {
            for my $ij (keys %{$file_keep_reads[$f]{$keep_read_id}}) {
                my ($i,$j) = split '_',$ij;
                print {$output_files[$i][$j]} $_;
            }
        }
    }
    close($fh);
    $f++;
}





__END__
