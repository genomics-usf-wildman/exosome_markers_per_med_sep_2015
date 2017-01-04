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
use IO::Uncompress::AnyUncompress;
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
if (not @ARGV) {
    push @USAGE_ERRORS,"You must provide at least one fasta file";
}
if ($options{paired} and (@files % 2 != 0)) {
    push @USAGE_ERRORS,"If reads are paired, you must give an even number of fasta files";
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
OUTER: for my $i ($options{paired}?(0..floor(@files/2)):(0..$#files)) {
    if ($options{paired}) {
        $i *= 2;
    }
    my $fh = open_fastq($files[$i]);
    my $fh2 = open_fastq($files[$i+1]) if $options{paired};
    while (not $fh->eof) {
        # are we likely to even keep this read?
        my $ent_match_rand = int rand($total_counts_matched+1);
        my $ent_unmatch_rand = int rand($total_counts_unmatched+1);
        if ($ent_match_rand > $max_samples and
            $ent_unmatch_rand > $max_samples
           ) {
            skip_a_read($fh);
            if ($options{paired}) {
                skip_a_read($fh2);
            }
            next;
        }
        my $read = read_a_read($fh);
        my $matched = $read->{read} =~ $motif_regex;
        my @entry = $read->{fastq};
        if ($options{paired}) {
            my $pair = read_a_read($fh2);
            if (not $matched) {
                $matched = $pair->{read} =~ $motif_regex;
            }
            push @entry,$pair->{fastq};
        }
        if ($matched) {
            $total_counts_matched++;
            reservoir_sampler($matched_samples,
                              $max_samples,
                              \@entry,
                              $total_counts_matched,
                              $ent_match_rand,
                             );
        } else {
            $total_counts_unmatched++;
            reservoir_sampler($unmatched_samples,
                              $max_samples,
                              \@entry,
                              $total_counts_unmatched,
                              $ent_unmatch_rand,
                             );
        }
    }
}
print STDERR "There were ".($total_counts_unmatched+$total_counts_matched).
    " reads in ".@files." files\n";
print STDERR $total_counts_matched." contained the motif, ".
    $total_counts_unmatched." did not\n";

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
        for (1..$read_count) {
            my ($matched_entry) = pop @{$matched_samples};
            my ($entry) = pop @{$unmatched_samples};
            if (rand(1) <= $options{motif_output_abundance}) {
                $entry = $matched_entry;
            }
            print {$output_files[$i][$j][0]} $entry->[0];
            if ($options{paired}) {
                print {$output_files[$i][$j][1]} $entry->[1];
            }
        }
    }
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
    my $fh;
    for (0..3) {
        my $gline = <$fh>;
    }
}

# read a single read from a fastq file
sub read_a_read {
    my ($fh) = @_;

    # read four lines
    my @read;
    my $read = '';
    while (<$fh>) {
        $read .= $_;
        chomp;
        push @read, $_;
        last if @read==4;
    }
    return if @read == 0;
    die "Truncated fastq file" if @read != 4;
    my ($read_id) = $read[0] =~ /^(\S+)/o;
    my ($quality) = sum(map {ord($_) - 33} split '',$read[3])/length($read[3]);
    return {quality => $quality,
            id => $read_id,
            read => $read[1],
            fastq => $read
           };
}

sub open_fastq {
    my ($fn) = @_;
    my $fh;
    open($fh,'-|','gzip','-dc',$fn) or
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
            for ($options{paired}?('_1','_2'):('')) {
                push @fns,
                    $options{output_prefix}.
                    '_r'.$options{reads}[$i].
                    '_s'.($j+1).$_.'.fastq';
            }
            for my $fn (@fns) {
                push @{$output_files[$i][$j]},
                    IO::File->new($fn,'w') or
                        die "Unable to open $fn for writing: $!";
            }
        }
    }
    return ($max_samples,@output_files);
}



__END__
