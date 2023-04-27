#!/usr/bin/perl

############################################################
#
# script to run geneid with hints.
#
# Author: J. Gómez Garrido
############################################################

use strict;
use warnings;

my $chunkfile = $ARGV[0];
my $junctions = $ARGV[1];
my $geneid_parameters = $ARGV[2];
my $geneid_options = $ARGV[3];
my $path = $ARGV[4];

my %scaff;
open JUNCTIONS, "<", "$junctions";
while (<JUNCTIONS>){
    chomp;
    my @line = split /\s+/, $_;
    $scaff{$line[0]}++ if (!exists $scaff{$line[0]});
    $scaff{$line[0]} .= "$_\n";
}
close JUNCTIONS;

open FASTA, "<", "$chunkfile";
my $id; 
my %keep; 
while (<FASTA>) {
    chomp;
    if ( $_ =~ m/\>/ ){
      my @ids = split '>', $_;
      $id = $ids[1];
    #  close SEQFILE if (SEQFILE);
      open HINTFILE, ">", "$id.hints.gff";
      print HINTFILE "$scaff{$id}";
      close HINTFILE;
      open SEQFILE, ">", "$id.masked.fa";
      $keep{$id}++;
    }
    print SEQFILE "$_\n";  
}
close FASTA;

foreach my $seq (keys %keep){
  print "Running: $path/bin/geneid -R $seq.hints.gff -P $geneid_parameters $geneid_options $seq.masked.fa  > $seq.geneid_introns.gff3\n";
  system "$path/bin/geneid -R $seq.hints.gff -P $geneid_parameters $geneid_options $seq.masked.fa  > $seq.geneid_introns.gff3";
}