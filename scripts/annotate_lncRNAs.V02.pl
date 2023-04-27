#!/usr/bin/perl

############################################################
#
# script to annotate the lncRNAs out of the pasa assemblies which have not been annotated as protein coding.
#
# Author: J. Gómez Garrido
############################################################

use strict;
use warnings;

my $annotation = $ARGV[0];
defined($annotation) || die ("##ERROR## This script requires an annotation output file as first argument\n");
my $pasa_assemblies = $ARGV[1];
defined($pasa_assemblies) || die ("##ERROR## This script requires a file with the pasa assemblies as second argument.\n");
my $project = $ARGV[2];
defined($project) || die ("##ERROR## This script requires a project name as third argument\n");
my $version = $ARGV[3];
defined($version) || die ("##ERROR## This script requires a version letter as last argument\n");

open ANNOTATION, "<", "$annotation";
my %annotated;
while (<ANNOTATION>) {
    chomp;
    next if /^\#/o;
    next if /^$/o;
    #print "$_\n";
    my @line = split /\s+/, $_;  

    if ($line[2] eq 'gene'){
        my $gene;
        if ($line[8] =~ m/ID=([^;]+)/) {$gene = $1;}
        $annotated{$line[0]}{$gene}->{start} = $line[3];
        $annotated{$line[0]}{$gene}->{end} = $line[4];
    }
}
close ANNOTATION;

open PASA, "<", "$pasa_assemblies";
my %assemblies;
my %assembl_exons;
my $i;
while (<PASA>) {
    chomp;
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\s+/, $_; 
    my $align;
    my $asmbl;
    if ($line[8] =~ m/ID=([^;]+)/) {$align = $1;}
    if ($line[8] =~ m/Target=([^\s]+)/) {$asmbl = $1;}
    $assemblies{$line[0]}{$asmbl} -> {start} = $line[3] if (!exists $assemblies{$line[0]}{$asmbl});
    $assemblies{$line[0]}{$asmbl} -> {end} = $line[4];
    if (!exists $assembl_exons{$asmbl}){
        $i = 1;
    }
    $assembl_exons{$asmbl}{$i} = $_;
    $i++;
}
close PASA;

my %coding;
foreach my $scaffold (sort keys %assemblies) {
    foreach my $pasa_asmbl (keys %{$assemblies{$scaffold}}) {
        my $start = $assemblies{$scaffold}{$pasa_asmbl}->{start};
        my $end = $assemblies{$scaffold}{$pasa_asmbl}->{end};
        foreach my $gene (keys %{$annotated{$scaffold}}) {
            if ($annotated{$scaffold}{$gene}->{start}<= $start && $annotated{$scaffold}{$gene}->{end} > $start) {
                #print "$scaffold\t$gene\t$annotated{$scaffold}{$gene}->{start}\t$annotated{$scaffold}{$gene}->{end}\t$pasa_asmbl\t$start\t$end\n";
                $coding{$pasa_asmbl}++ if (!exists $coding{$pasa_asmbl});
            }
            elsif ($annotated{$scaffold}{$gene}->{start} > $start && $annotated{$scaffold}{$gene}->{end} <= $end) {
                #print "$scaffold\t$gene\t$annotated{$scaffold}{$gene}->{start}\t$annotated{$scaffold}{$gene}->{end}\t$pasa_asmbl\t$start\t$end\n";
                $coding{$pasa_asmbl}++ if (!exists $coding{$pasa_asmbl});
            }
            elsif ($annotated{$scaffold}{$gene}->{start} < $end && $annotated{$scaffold}{$gene}->{end} > $end) {
               # print "$scaffold\t$gene\t$annotated{$scaffold}{$gene}->{start}\t$annotated{$scaffold}{$gene}->{end}\t$pasa_asmbl\t$start\t$end\n";
                $coding{$pasa_asmbl}++ if (!exists $coding{$pasa_asmbl});
            }            
        }
    }
}

$i = 1;
foreach my $scaffold (sort keys %assemblies) {
    foreach my $pasa_asmbl (sort keys %{$assemblies{$scaffold}}) {
        my $length = $assemblies{$scaffold}{$pasa_asmbl}->{end} - $assemblies{$scaffold}{$pasa_asmbl}->{start};
        if (!exists $coding{$pasa_asmbl} && $length >= "200") {            
            my $lncrna_id = $project . "lnc" . $version . sprintf("%06d",$i);
            my @line = split /\s+/, $assembl_exons{$pasa_asmbl}{1};
            print "$scaffold\tCNAG\ttranscript\t$assemblies{$scaffold}{$pasa_asmbl}->{start}\t$assemblies{$scaffold}{$pasa_asmbl}->{end}\t.\t$line[6]\t.\tID=$lncrna_id;Name=$lncrna_id;Target=$pasa_asmbl;Description=lncRNA\n";
            foreach my $exon (sort {$a<=>$b} keys %{$assembl_exons{$pasa_asmbl}}) {
                my @line = split /\s+/, $assembl_exons{$pasa_asmbl}{$exon};
                print "$scaffold\tCNAG\texon\t$line[3]\t$line[4]\t.\t$line[6]\t.\tParent=$lncrna_id;ID=$lncrna_id" . ".exon" . "$exon;Name=$lncrna_id" . ".exon\n";
                print "$scaffold\tCNAG\tCDS\t$line[3]\t$line[4]\t.\t$line[6]\t.\tParent=$lncrna_id;ID=$lncrna_id" . "_cds;Name=$lncrna_id" . "_cds\n";
            }
            $i++;
        }
      #  print "$assembl_exons{$pasa_asmbl}" if (!exists $coding{$pasa_asmbl});
    }
}
