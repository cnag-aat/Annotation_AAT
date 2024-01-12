#!/usr/bin/perl

############################################################
#
# script to select the best EVM output after running it with several weights.
#
# Author: J. Gómez Garrido
############################################################
use strict;
use File::Basename;

my $reference_transcripts = $ARGV[0];

my %bp;

for (my $i=1; $i < scalar @ARGV; $i++) {
   my $name = basename($ARGV[$i], ".gff3");
   system "gawk \'\$3==\"CDS\"\' $ARGV[$i] > $name.CDS.gff3";
   system "intersectBed -wo -a $name.CDS.gff3 -b $reference_transcripts  > $name.BT.gmap.out";
   open BEDTOOLS, "<", "$name.BT.gmap.out" ;
   $bp{$i} = 0;
   while (<BEDTOOLS>) {
     chomp;
     my @line = split /\t/, $_;
     my $current = pop @line;
     $bp{$name} = $bp{$i} + $current;
   }
   close BEDTOOLS; 
}

my $best = 0;
my $bestid = "";
foreach (keys %bp) {
   if ($bp{$_} > $best) {
      $best = $bp{$_};
      $bestid = $_;
   }
}

system "ln -s $bestid.gff3 evm.best.gff3";
print "Best evm weights are $bestid\n";
