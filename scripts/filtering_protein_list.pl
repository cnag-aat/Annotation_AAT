#!/usr/bin/env perl
use strict;
#use warnings;
use Getopt::Long;

my ($gfffile, $ids_list,$add, $discard);

GetOptions(
           'gff:s'           => \$gfffile, 
           'l:s'                => \$ids_list,
           'a'              => \$add,   #if this option is given the lines with the id are going to be kept
           'b'               => \$discard, #If this option is given, the lines with the id are going to be discarded
           );
           
open IDS,"< $ids_list" || die "cannot open $ids_list";
my %ids;
my %genes;
while(<IDS>){
  chomp;
  $ids{$_}++;
 # print "$gene\n";
}
close IDS;

open GFF, "<", "$gfffile" || die "Cannot open $gfffile";

my $protein;
my %seen;
my %trans;
my $tr;
if ($add) {
  while (<GFF>) {
    chomp;
    next if /^$/o;
    if ($_ =~ m/#/){
      print "$_\n";
      next;
    }
    my @line = split /\t/, $_;  
    if ($line[2] eq 'gene') {
      my $gene;
      if ($line[8] =~ m/ID=([^;]+)/) {
        $gene = $1;
      }
      elsif ($line[8] =~ m/ID=/) {
        $gene = $';
      }
      $genes{$gene} = $_;
    }
    elsif ($line[2] =~ m/transcript|mRNA/) {
      if ($line[8] =~ m/product=([^;]+)/) {
        $protein = $1;
      }
      elsif ($line[8] =~ m/product=/) {
        $protein = $';
      }  
      if ($line[8] =~ m/ID=([^;]+)/) {$tr = $1;}
      elsif ($line[8] =~ m/ID=/) {$tr = $';}  
      if (exists $ids{$protein}){
        my $gene = $tr;
        $gene =~ s/T([0-9]*)//;
      #  print "$gene\n";
        print "$genes{$gene}\n" if (!exists $seen{$gene});
        $seen{$gene}++;
        print "$_\n";
        $trans{$tr}++;
      }
    }
    else {
      if ($line[8] =~ m/Parent=([^;]+)/) {$tr = $1;}
      elsif ($line[8] =~ m/Parent=/) {$tr = $';}
      print "$_\n" if (exists $trans{$tr});
    }
  }
}

if ($discard) {
  while (<GFF>) {
    chomp;
    next if /^$/o;
    if ($_ =~ m/#/){
      print "$_\n";
      next;
    }
    my @line = split /\t/, $_;  
    if ($line[2] eq 'gene') {
      my $gene;
      if ($line[8] =~ m/ID=([^;]+)/) {
        $gene = $1;
      }
      elsif ($line[8] =~ m/ID=/) {
        $gene = $';
      }
      $genes{$gene} = $_;
    }
    elsif ($line[2] =~ m/transcript|mRNA/) {
      if ($line[8] =~ m/product=([^;]+)/) {
        $protein = $1;
      }
      elsif ($line[8] =~ m/product=/) {
        $protein = $';
      }  
      if ($line[8] =~ m/ID=([^;]+)/) {$tr = $1;}
      elsif ($line[8] =~ m/ID=/) {$tr = $';}  
      if (!exists $ids{$protein}){
        my $gene = $tr;
        $gene =~ s/T([0-9]*)//;
      #  print "$gene\n";
        print "$genes{$gene}\n" if (!exists $seen{$gene});
        $seen{$gene}++;
        print "$_\n";
        $trans{$tr}++;
      }
    }
    elsif ($line[2] =~ m/exon|CDS/){
      if ($line[8] =~ m/Parent=([^;]+)/) {$tr = $1;}
      elsif ($line[8] =~ m/Parent=/) {$tr = $';}
      print "$_\n" if (exists $trans{$tr});
    }
  }
}
