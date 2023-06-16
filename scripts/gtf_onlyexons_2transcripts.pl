#!/usr/bin/env perl

use strict;
use warnings;

my %trans; 
my %exons;
while (<STDIN>){
  chomp;
  my @line = split /\t/, $_;
  next if ($line[2] ne "exon");
  my ($gene, $trans);
  my $exon=1;
  if ($_ =~ m/gene_id \"([^\"]+)/){
    $gene = $1;
  }
  if ($_=~ m/transcript_id \"([^\"]+)/){
    $trans = $1;
  }
  if (!exists $trans{$trans}){
    $trans{$trans}->{start} = $line[3];
    $trans{$trans}->{end} = $line[4];
  }
  else {
    $trans{$trans}->{start}=$line[3] if ($trans{$trans}->{start} > $line[3]);
    $trans{$trans}->{end}=$line[4] if ($trans{$trans} ->{end} < $line[4]);
    while (exists $exons{$gene}{$trans}{$exon}){
      $exon++;
    }
  }
  if ($_ =~m/exon_number \"([^\"]+)/){
    $exon = $1;
  }
  else {
    $_ .= "exon_number \"$exon\";";
  }
  $exons{$gene}{$trans}{$exon}=$_;
}

my %seen;
foreach my $g (sort keys %exons){
  foreach my $t (sort keys %{$exons{$g}}){
    foreach my $e (sort {$a<=>$b} keys %{$exons{$g}{$t}}){
      my $line = $exons{$g}{$t}{$e};
      if ($e == 1){
        my @line = split /\t/, $line;
        my @trans_line = @line;
        $trans_line[2] = "transcript";
        $trans_line[3] = $trans{$t}->{start};
        $trans_line[4] = $trans{$t}->{end};
        $trans_line[8] =~ s/exon_number \"$e\"; //;
        my $trans_line = join "\t", @trans_line;
        print "$trans_line\n";
      }
      print "$line\n";
    }
  }
}
