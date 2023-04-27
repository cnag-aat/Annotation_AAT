#!/usr/bin/perl

############################################################
#
# script to select the best EVM output after running it with several weights.
#
# Author: J. Gómez Garrido
############################################################
use strict;

my $reference_transcripts = $ARGV[0];
my $total_evm = $ARGV[1];
#print "$gmap_alignments\n";

my %bp;
my $i;
my @line;

for ($i=1; $i<=$total_evm; $i++) {
   #print "$i\n";
   system "gawk \'\$3==\"CDS\"\' evm_weights_'$i'.gff3 > evm_weights_'$i'.CDS.gff3";
   system "intersectBed -wo -a evm_weights_'$i'.CDS.gff3 -b $reference_transcripts  > evm_weights_'$i'.BT.gmap.out";
   open BEDTOOLS, "<", "evm_weights_$i.BT.gmap.out" ;
   $bp{$i} = 0;
   while (<BEDTOOLS>) {
      chomp;
      @line = split /\t/, $_;
      my $current = pop @line;
     # print "$current\n";
      $bp{$i} = $bp{$i} + $current;
    #  print "$bp{$i}\n";
   }
  # print "$bp{$i}\n";
   close BEDTOOLS; 
 
}

my $best = 0;
my $bestid = "";
foreach (keys %bp) {
  # print "$_\n";
  # print "$bp{$_}\n";
  # print "$best\n";
   if ($bp{$_} > $best) {
      $best = $bp{$_};
      $bestid = $_;
   }
}

system "ln -s evm_weights_'$bestid'.gff3 evm.best.gff3";
print "Best evm weights are $bestid\n";
