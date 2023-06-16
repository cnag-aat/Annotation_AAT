#!/usr/bin/env perl

use strict;
use warnings;

my %trans; 
my %exons;
my $exon;
while (<STDIN>){
  chomp;
  my @line = split /\t/, $_;
  next if ($line[2] ne "cDNA_match");
  my ($gene, $trans);
  if ($_ =~ m/ID=([^;]+)/){
    $gene = $1;
    $trans = $1 . ".mrna";
  }
  if (!exists $trans{$trans}){
    $trans{$trans}->{start} = $line[3];
    $trans{$trans}->{end} = $line[4];
    $exon = 1;
  }
  else {
    $trans{$trans}->{start}=$line[3] if ($trans{$trans}->{start} > $line[3]);
    $trans{$trans}->{end}=$line[4] if ($trans{$trans} ->{end} < $line[4]);
  }
  $exons{$gene}{$trans}{$exon}=$_;
  $exon++;
}

my %seen;
foreach my $g (sort keys %exons){
  foreach my $t (sort keys %{$exons{$g}}){
    foreach my $e (sort {$a<=>$b} keys %{$exons{$g}{$t}}){
      my $line = $exons{$g}{$t}{$e};
      my @line = split /\t/, $line;
      if ($e == 1){
        my @trans_line = @line;
        $trans_line[2] = "mRNA";
        $trans_line[3] = $trans{$t}->{start};
        $trans_line[4] = $trans{$t}->{end};
        $trans_line[8] = "ID=$t;";
        my $trans_line = join "\t", @trans_line;
        print "$trans_line\n";
      }
      $line[2] = "exon";
      $line[8] = "ID=$t.$e;Parent=$t";
      $line = join "\t", @line;
      print "$line\n";
    }
  }
}

