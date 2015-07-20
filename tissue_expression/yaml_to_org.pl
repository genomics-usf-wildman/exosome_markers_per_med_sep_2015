#!/usr/bin/perl
# yaml_to_org.pl coverts a yaml file into an org table
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2014 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

yaml_to_org.pl - coverts a yaml file into an org table

=head1 SYNOPSIS

yaml_to_org.pl [options]

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

yaml_to_org.pl

=cut


use vars qw($DEBUG);

use YAML;

my %options = (debug           => 0,
               help            => 0,
               man             => 0,
              );

GetOptions(\%options,
           'include=s@',
           'exclude=s@',
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;

if (@ARGV != 1) {
    push @USAGE_ERRORS,'You must provide exactly one yaml file';
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;

my $yaml = YAML::LoadFile($ARGV[0]);

my @columns;

if (not defined $options{include}) {
    my %cols = map {map {($_,1)} keys %{$_}} @{$yaml};
    @columns = keys %cols;
} else {
   @columns = @{$options{include}};
}

if (defined $options{exclude}) {
    my %exc = map {($_,1)} @{$options{exclude}};
    @columns = grep {not exists $exc{$_}} @columns;
}

binmode(STDOUT,":utf8");

print "|--|--|\n";
print '|'.join('|',@columns)."|\n";
print "|--|--|\n";
for my $y (@{$yaml}) {
    print '|'.join('|',map{defined $y->{$_}?$y->{$_}:''} @columns)."|\n";
}
print "|--|--|\n";


__END__
