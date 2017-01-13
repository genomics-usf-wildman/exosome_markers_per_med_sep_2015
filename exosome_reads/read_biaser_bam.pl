#!/usr/bin/env perl
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
   --paired Are reads paired?
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
use vars qw($DEBUG);
#use Math::Random;
use List::Util qw(shuffle);
use POSIX qw(floor);

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
               paired          => 0,
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
           'paired!',
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
# random_set_seed(random_seed_from_phrase($options{seed}));
print STDERR "Seed: $options{seed}\n";

my @files = @ARGV;
if (@ARGV != 1) {
    push @USAGE_ERRORS,"You must provide exactly one bam file";
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;

# build motif markers regex
my $motif_regex = '('.join('|',map {s/[UT]/\[TU\]/; $_} @{$options{motif}}).')';
print STDERR "motif regex: $motif_regex\n";
$motif_regex = qr/$motif_regex/;

my ($max_samples,@output_files) = open_output_files(%options);

## read the fastq files and store randomly selected reads in the
## reservoir sampler
my $total_counts_matched = 0;
my $total_counts_unmatched = 0;
my $matched_samples = [];
my $unmatched_samples = [];
my $header = '';
my $in_header = 1;
my $fh = open_bam($files[0]);
while (<$fh>) {
    if ($in_header and /^@/) {
        $header .= $_;
        next;
    } else {
        $in_header = 0;
    }
    # are we likely to even keep this read?
    my $ent_match_rand = int rand($total_counts_matched+1);
    my $ent_unmatch_rand = int rand($total_counts_unmatched+1);
    if ($ent_match_rand > $max_samples and
        $ent_unmatch_rand > $max_samples
       ) {
        next;
    }
    my @r = split /\t/;
    my $matched = $r[9] =~ $motif_regex;
    if ($matched) {
        $total_counts_matched++;
        reservoir_sampler($matched_samples,
                          $max_samples,
                          \@r,
                          $total_counts_matched,
                          $ent_match_rand,
                         );
    } else {
        $total_counts_unmatched++;
        reservoir_sampler($unmatched_samples,
                          $max_samples,
                          \@r,
                          $total_counts_unmatched,
                          $ent_unmatch_rand,
                         );
    }
}
print STDERR "There were ".($total_counts_unmatched+$total_counts_matched).
    " reads in ".@files." files\n";
print STDERR $total_counts_matched." contained the motif, ".
    $total_counts_unmatched." did not\n";
close($fh) or die "Unable to finish reading bam: $!";

# shuffle the matching and unmatched reads
@{$matched_samples} = shuffle(@{$matched_samples});
while (@{$matched_samples} < $max_samples) {
    push @{$matched_samples},@{$matched_samples};
}
@{$unmatched_samples} = shuffle(@{$unmatched_samples});
while (@{$unmatched_samples} < $max_samples) {
    push @{$unmatched_samples},@{$unmatched_samples};
}

# output sampled reads into the output files
for my $i (0..$#output_files) {
    my $read_count = $options{reads}[$i];
    for my $j (0..$#{$output_files[$i]}) {
        my @outputs;
        for (1..$read_count) {
            my ($matched_entry) = pop @{$matched_samples};
            my ($entry) = pop @{$unmatched_samples};
            if (rand(1) <= $options{motif_output_abundance}) {
                $entry = $matched_entry;
            }
            push @outputs,$entry;
        }
        # open output file
        output_sampled_reads(\@outputs,$header,$output_files[$i][$j]);
    }
}

sub output_sampled_reads {
    my ($outputs,$header,$fn) = @_;

    my $fh;
    open($fh,'|-','samtools','view','-b','-o',$fn,'-') or
        die "Unable to open samtools view -b -o $fn - for writing: $!";
    print {$fh} $header;
    for my $entry (sort { $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3] }
                   @{$outputs}
                  ) {
        print {$fh} join("\t",@{$entry});
    }
    close($fh) or die "Unable to close samtools view: $!";
}

# Ideally we would use a weighted random sampler, but because we only
# know the output abundance and not the input abundance, we cannot use
# weights
sub reservoir_sampler {
    my ($samples,$max_samples,$entry,$n_entries,$ent_rand) = @_;
    if (@{$samples} < $max_samples) {
        push @{$samples},$entry;
        return;
    }
    if ($ent_rand <= $max_samples) {
        $samples->[floor(rand(@{$samples}))] =
            $entry;
    }
}

sub skip_a_read {
    my $fh = (@_);
    my $gline = <$fh>;
}

sub open_bam {
    my ($fn) = @_;
    my $fh;
    open($fh,'-|','samtools','view','-h',$fn) or
        die "Unable to open $fn for reading: $!";
    return $fh;
}

sub open_output_files {
    my %options = @_;
    my $max_samples = 0;
    my @output_files;
    for my $i (0..$#{$options{reads}}) {
        for my $j (0..($options{samplings}-1)) {
            $max_samples += $options{reads}[$i];
            my @fns;
            $output_files[$i][$j] =
                $options{output_prefix}.
                '_r'.$options{reads}[$i].
                '_s'.($j+1).'.bam';
        }
    }
    return ($max_samples,@output_files);
}



__END__
