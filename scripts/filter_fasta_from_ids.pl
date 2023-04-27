#! /usr/bin/perl

use warnings;
use strict;

use Getopt::Long;

my ($fasta, $ids, $add, $discard);

GetOptions(
           'f:s'           => \$fasta, 
           'l:s'                => \$ids,
           'a'              => \$add,   #if this option is given the lines with the id are going to be kept
           'b'               => \$discard, #If this option is given, the lines with the id are goin to be discarded
           );


open IDS, "<", "$ids";
my %ids;
while (<IDS>){
  chomp;
  my @line = split /\s+/, $_;
  $ids{$line[0]}++;
}
close IDS;

open FASTA, "<", "$fasta";
my $seq;
while (<FASTA>){
  chomp;
  next if ($_ !~ m/./);
 # if ($_ =~ m/>([^|]+)/ || $_ =~ m/>([^l]+)/){
   if ($_ =~ m/>([^ ]+)/){
    $seq = $1;
    #print "$seq\n";
  }
  if ($add){
    print "$_\n" if (exists $ids{$seq});
  }
  else {
    print "$_\n" if (!exists $ids{$seq});
  }
}
close FASTA;

