#!/usr/bin/perl

use warnings;
use strict;

use utf8;

use XML::LibXML;
use Data::Printer;

my $dom = XML::LibXML->load_xml(location => $ARGV[0]);

# print map {$_->nodeName."\n"} map {$_->childNodes() } $dom->childNodes();

my @samples = $dom->firstChild->getChildrenByLocalName('Sample');

#p @samples;

# sub get_content 

sub fixup_text{
    my ($text) = @_;
    $text =~ s/[\t\n\s]/ /g;
    return($text);
}
binmode(STDOUT,':utf8');

print "title\taccession\tgrowth\tsra\n";
for my $sample (@samples) {
    # print $sample->nodeName()."\n";
    for (qw(Title Accession Growth-Protocol)) {
        my @t = $sample->findnodes(".//*[local-name()='$_']");
        if (not defined $t[0]) {
            print "\t";
        } else {
            print fixup_text($t[0]->textContent())."\t";
        }
    }
    my @t = $sample->findnodes(".//*[local-name()='Supplementary-Data']");
    for my $t (@t) {
        # p $t;
        if (0 < scalar grep {$_->nodeName eq 'type' and
                                 $_->value eq 'SRA Experiment'}
            $t->attributes()) {
            print fixup_text($t->textContent());
        }
    }
    print "\n";
}


