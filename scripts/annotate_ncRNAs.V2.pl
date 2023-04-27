#!/usr/bin/perl

############################################################
#
# script to get a small ncRNA annotation taking the outputs of infernal and tRNAscan-SE.
#
# Author: J. Gómez Garrido
############################################################

use strict;
use warnings;

my $infernal = $ARGV[0];
defined($infernal) || die ("##ERROR## This script requires an infernal cmsearch tabular output as argument\n");

my $trnascan = $ARGV[1];
defined($trnascan) || die ("##ERROR## This script requires a tRNAscan-SE output as second argument\n");

my $project = $ARGV[2];
defined($project) || die ("##ERROR## This script requieres a project name as third argument\n");

my $version = $ARGV[3];
defined($version) || die ("##ERROR## This script requieres a version letter as last argument\n");

my %scaff;
my %exons;

open INFERNAL, "<", "$infernal";
while (<INFERNAL>) {
    next if /^\#/o;
    chomp;
    my @line = split /\s+/, $_;
   # print "$line[17]\n";
    my $scaffold = $line[0];
    my $rnatype = $line[2];    
    my $modelname = $line[3];
    my $seqstart = $line[7];
    my $seqend = $line[8];
    my $strand = $line[9];
    my $trunc = $line[10];
    my $score = $line[14];
    my $eval = $line[15];
    if ($rnatype !~ m/tRNA/) {
        if ($eval < 0.001 && $strand eq '+') {
            $scaff{$scaffold}{$seqstart}= "$scaffold\tCNAG\tncRNA\t$seqstart\t$seqend\t$score\t$strand\t.\tID=" . $project . "snc" . $version . ";Type=$rnatype;RFAM_model=$modelname\n";
            $exons{$scaffold}{$seqstart} = "$scaffold\tCNAG\texon\t$seqstart\t$seqend\t.\t$strand\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon1\n";
        }
        elsif ($eval < 0.001 && $strand eq '-') {
            $scaff{$scaffold}{$seqend}= "$scaffold\tCNAG\tncRNA\t$seqend\t$seqstart\t$score\t$strand\t.\tID=" . $project . "snc" . $version . ";Type=$rnatype;RFAM_model=$modelname\n";
            $exons{$scaffold}{$seqend} = "$scaffold\tCNAG\texon\t$seqend\t$seqstart\t.\t$strand\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon1\n";
        }
      #  print "$scaffold\tINFERNAL\tncRNA\t$seqstart\t$seqend\t$score\t$strand\t.\tID=;Type=$rnatype;RFAM_model=$modelname\n" if ($eval < 0.001 && $strand eq '+');
       # print "$scaffold\tINFERNAL\tncRNA\t$seqend\t$seqstart\t$score\t$strand\t.\tID=;Type=$rnatype;RFAM_model=$modelname\n" if ($eval < 0.001 && $strand eq '-');
    }
}
close INFERNAL;

open tRNA, "<", "$trnascan";
while (<tRNA>) {
    chomp;
    next if /^\#/o;
    my @line = split /\s+/, $_;
   # print "$line[17]\n";
    my $scaffold = $line[0];
    my $start = $line[2];
    my $end = $line[3];
    my $intronstart = $line[6];
    my $intronend = $line[7];
    my $score = $line[8];
    my $trnatype = $line[4];
    next if ($trnatype eq "Pseudo");
    if ($start < $end) {
       # print "$scaffold\ttRNAscan-SE\tncRNA\t$start\t$end\t$score\t+\t.\tID=;Type=tRNA.$trnatype;\n";
         $scaff{$scaffold}{$start}= "$scaffold\tCNAG\tncRNA\t$start\t$end\t$score\t+\t.\tID=" . $project . "snc" . $version . ";Type=tRNA.$trnatype\n";
         if ($intronstart > 0) {
             my $exonend = $intronstart - 1;
             my $exonstart = $intronend + 1;
             $exons{$scaffold}{$start} = "$scaffold\tCNAG\texon\t$start\t$exonend\t.\t+\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon1\n$scaffold\tCNAG\texon\t$exonstart\t$end\t.\t+\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon2\n";
         }
         else {
             $exons{$scaffold}{$start} = "$scaffold\tCNAG\texon\t$start\t$end\t.\t+\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon1\n";
         }
    }
    else {
       # print "$scaffold\ttRNAscan-SE\tncRNA\t$end\t$start\t$score\t-\t.\tID=;Type=tRNA.$trnatype;\n";
         $scaff{$scaffold}{$end}= "$scaffold\tCNAG\tncRNA\t$end\t$start\t$score\t-\t.\tID=" . $project . "snc" . $version . ";Type=tRNA.$trnatype\n";
         if ($intronend > 0) {
             my $exonend = $intronend - 1;
             my $exonstart = $intronstart + 1;
             $exons{$scaffold}{$end} = "$scaffold\tCNAG\texon\t$end\t$intronend\t.\t-\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon2\n$scaffold\tCNAG\texon\t$intronstart\t$start\t.\t-\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon1\n";
         }
         else {
             $exons{$scaffold}{$end} = "$scaffold\tCNAG\texon\t$end\t$start\t.\t-\t.\tParent=" . $project . "snc" . $version . ";ID=" . $project . "snc" . $version . ".exon1\n";
         }
    }
}
close tRNA;

my $i = 1;
foreach my $seq (sort keys %scaff) {
    foreach (sort keys %{$scaff{$seq}}) {
        my $rna_id = sprintf("%06d",$i);
        $scaff{$seq}{$_} =~ s/;Type=/$rna_id;Type=/;
        $exons{$seq}{$_} =~ s/;ID=/$rna_id;ID=/g if (exists $exons{$seq}{$_});
        $exons{$seq}{$_} =~ s/\.exon/$rna_id.exon/g if (exists $exons{$seq}{$_});
        print "$scaff{$seq}{$_}";
        print "$exons{$seq}{$_}" if (exists $exons{$seq}{$_});
        $i++;
        
    } 
}
