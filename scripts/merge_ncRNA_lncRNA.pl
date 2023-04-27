#!/usr/bin/perl

############################################################
#
# script to get a consensus ncRNA annotation merging the small ncRNA annotation and the lncRNA annotation.
#
# Author: J. Gómez Garrido
############################################################

use strict;
use warnings;

my $ncRNA = $ARGV[0];
defined($ncRNA) || die ("##ERROR## This script requires a ncRNA annotation file as first argument\n");

my $lncRNA = $ARGV[1];
defined($lncRNA) || die ("##ERROR## This script requires a lncRNA annotation file as second argument\n");

my $project = $ARGV[2];
defined($project) || die ("##ERROR## This script requires a project name with the assembly version as third argument\n");

my $nc_version = $ARGV[3];
defined($nc_version) || die ("##ERROR## This script requires a version of the nc annotation as last argument\n");


open ncRNA, "<", "$ncRNA";
my %ncRNA;
my %ncRNA_exons;
while (<ncRNA>){
    chomp;
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\t/, $_;
    if ($line[2] eq 'ncRNA') {
        my $id;
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        $ncRNA{$line[0]}{$line[3]}{$line[6]}{nc}{$id} = $_;
      #  print "$line[6]\n";
    }
    elsif ($line[2] eq 'exon') {
        my $parent;
        my $id;
        $line[8] .= ";";
        if ($line[8] =~ m/Parent=([^;]+)/) {$parent = $1;}
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        $ncRNA_exons{$parent}{$id} = $_;
    }
}
close ncRNA;

open lncRNA, "<", "$lncRNA";
my %genes;
my %throw;
my %lnctranscript;
my $gene;
while (<lncRNA>) {
    chomp;
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\t/, $_;   
    if ($line[2] eq 'gene') {
        if ($line[8] =~ m/ID=/) {$gene = $';}
        $ncRNA{$line[0]}{$line[3]}{$line[6]}{lnc}{$gene} = $_;        
    }
    elsif ($line[2] =~ m/ncRNA|transcript/) {
        my $id;
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        my @tname = split /T/, $id;
        pop @tname;
        $gene = join "T", @tname;
       # print "!!$gene\t$id\t@tname\n";
        #my @gene = split /T/, $id;
        if ($line[8] =~ m/small_ncRNA/) {$throw{$gene}++;}
        $lnctranscript{$gene}{$id} = $_;
    }
    elsif ($line[2] eq 'exon') {
        my $parent;
        my $id;
        if ($line[8] =~ m/Parent=([^;]+)/) {$parent = $1;}
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        $ncRNA_exons{$parent}{$id} = $_;
    }
}
close lncRNA;

print "##gff-version 3\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "# date: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

my $i = 1;
foreach my $scaffold (sort keys %ncRNA) {
    my $gene_id;
    foreach my $start (sort {$a<=>$b} keys %{$ncRNA{$scaffold}}) {
      foreach my $strand (keys %{$ncRNA{$scaffold}{$start}}){
       # print "!!$scaffold\t$start\t$strand\n";
        if (exists $ncRNA{$scaffold}{$start}{$strand}{nc}) { 
            foreach my $transcript (keys %{$ncRNA{$scaffold}{$start}{$strand}{nc}}) {
                $gene_id = $project . "nc" . $nc_version . sprintf("%06d",$i);
                my @line = split /\t/, $ncRNA{$scaffold}{$start}{$strand}{nc}{$transcript};
                $line[8] .= ";";
                my $type;
                if ($line[8] =~ m/Type=([^;]+)/ || $line[8] =~ m/TYPE=([^;]+)/ ) {$type = $1;}
                print "\#\#\#\n$scaffold\tCNAG\tgene\t$start\t$line[4]\t.\t$line[6]\t$line[7]\tID=$gene_id\n";
                print "$scaffold\tCNAG\ttranscript\t$start\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tID=$gene_id" . "T1;Parent=$gene_id;Name=$gene_id" . "T1;Description=$type\n";
                foreach my $exons (sort keys %{$ncRNA_exons{$transcript}}) {
                    my @line = split /\t/, $ncRNA_exons{$transcript}{$exons};
                    my $exon;
                    $line[8] = $line[8] . ";";
                    if ($line[8] =~ m/exon([^;]+)/) {$exon = $1;}
                    print "$scaffold\tCNAG\texon\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tParent=$gene_id" . "T1;ID=$gene_id" . "T1.exon" . "$exon;Name=$gene_id" . "T1\n";
                }
            }
            $i++;
        }
        if (exists $ncRNA{$scaffold}{$start}{$strand}{lnc}) {
            foreach my $gene (keys %{$ncRNA{$scaffold}{$start}{$strand}{lnc}}) {
              #  print "!!$gene\t$strand\n";
                if (!exists $throw{$gene}) {
                    $gene_id = $project . "nc" . $nc_version . sprintf("%06d",$i);  
                    my @line = split /\t/, $ncRNA{$scaffold}{$start}{$strand}{lnc}{$gene};
                    $line[8] = "ID=$gene_id";
                    local $" = "\t";
                    print "\#\#\#\n@line\n";  
                    foreach my $transcript (sort keys %{$lnctranscript{$gene}}) { 
                        my $t = 1;
                        my $old_id;
                        my $attributes;
                        my @line = split /\t/, $lnctranscript{$gene}{$transcript};
                        $line[8] .= ";";
                        $line[8] =~ s/ID=([^;]+)//;
                        $line[8] =~ s/Parent=([^;]+)//;
                        if ($line[8] =~ m/Name=([^;]+)/) {$old_id = $1; $attributes = $'}                 
                        if ($old_id =~ m/T/) {
                          my @tname = split /T/, $old_id;
                          my $index = scalar @tname - 1;
                          $t = $tname[$index];
                        }                          
                        print "$scaffold\tCNAG\ttranscript\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tID=$gene_id" . "T" . "$t;Parent=$gene_id;Name=$gene_id" . "T" . "$t" . "$attributes\n";
                        foreach my $exons (sort keys %{$ncRNA_exons{$transcript}}) {
                            my @line = split /\t/, $ncRNA_exons{$transcript}{$exons};
                            my $exon;
                            if ($line[8] =~ m/\.exon([^;]+)/) {$exon = $1;}
                            print "$scaffold\tCNAG\texon\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tParent=$gene_id" . "T" . "$t;ID=$gene_id" . "T" . "$t.exon" . "$exon;Name=$gene_id" . "T" . "$t\n";
                        }
                    }
                }
                $i++;
            }
        }
      }
    }
}
