#!/usr/bin/env perl
use strict;
#use warnings;
use Getopt::Long;

my ($gfffile, $ids_list,$add, $discard);

GetOptions(
           'gff:s'           => \$gfffile, 
           'l:s'                => \$ids_list,  
           'a'              => \$add,   #if this option is given the lines with the id are going to be kept
           'b'               => \$discard, #If this option is given, the lines with the id are goin to be discarded
           );
           
open IDS,"< $ids_list" || die "cannot open $ids_list";
my %ids;
while(<IDS>){
  chomp;
  $ids{$_}++;
}
close IDS;

open GFF, "<", "$gfffile" || die "Cannot open $gfffile";

my @F8;

my $tr;
my %genes;
my %seen;
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
        if ($line[8] =~ m/ID=([^;]+)/) {$gene = $1;}
        elsif ($line[8] =~ m/ID=/) {$gene = $';}
        $genes{$gene} = $_;
        #print "$gene\n";
      }
      elsif ($line[2] =~ m/transcript|mRNA|ncRNA|tRNA|rRNA/) {
        if ($line[8] =~ m/ID=([^;]+)/) {$tr = $1;}
        elsif ($line[8] =~ m/ID=/) {$tr = $';} 
        my $gene;       
        if (exists $ids{$tr}){
          if ($_ =~ m/Parent=([^;]+)/){
            $gene = $1;
          }
          elsif ($_ =~m/Parent=/){
            $gene = $';
          }
          else {
            $gene = $tr;
            $gene =~ s/T([0-9]*)//;
          }
        #  print "$gene\n";
          print "$genes{$gene}\n" if (!exists $seen{$gene});
          $seen{$gene}++;
          print "$_\n";
        }
      }
    #  elsif ($line[2] =~ m/exon|CDS/){
      else{
        if ($line[8] =~ m/Parent=([^;]+)/) {$tr = $1;}
        elsif ($line[8] =~ m/Parent=/) {$tr = $';}
        print "$_\n" if (exists $ids{$tr});
      }
   }
}
my $gene;
if ($discard) {
  while (<GFF>) {
    chomp $_;
    if ($_ =~ m/#/){
      print "$_\n";
      next;
    }
    $_ = $_ . ";";
    my @line = split /\s+/, $_;  
    if ($line[2] eq 'gene') {
      if ($line[8] =~ m/ID=([^;]+)/) {$gene = $1;}
      elsif ($line[8] =~ m/ID=/) {$gene = $';}
     # print STDERR "!!repeated gene $gene\n";
      $genes{$gene} = $_;      
    } 
    elsif ($line[2] =~ m/transcript|mRNA|ncRNA|tRNA|rRNA/) {
      if ($line[8] =~ m/ID=([^;]+)/) {$tr = $1;}
      elsif ($line[8] =~ m/ID=/) {$tr = $';}        
      next if (exists $ids{$tr});
     # $gene = $tr;
     # $gene =~ s/T([0-9]*)//;
     # print "!!$gene\n" if (!exists $seen{$gene});
      print "$genes{$gene}\n" if (!exists $seen{$gene});
      $seen{$gene}++;
      print "$_\n";
    }
    elsif ($line[2] =~ m/exon|CDS/){
      if ($line[8] =~ m/Parent=([^;]+)/) {$tr = $1;}
      elsif ($line[8] =~ m/Parent=/) {$tr = $';}
      print "$_\n" if (!exists $ids{$tr});
    }
  }
}
close GFF;
