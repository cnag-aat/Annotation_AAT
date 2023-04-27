#!/usr/bin/env perl

use strict;
#use warnings;
use lib "/home/groups/assembly/talioto/myperlmods/";
#use lib "/apps/BIOPERL/1.6.1/lib/perl5";
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;
use File::Basename qw( fileparse );
use SeqOp;

my $gfffile =0; 
my $fa= 0;

GetOptions(
           'i|g|file|gff:s'             => \$gfffile,
           'f|s|seq:s'          => \$fa
          );

open (IN, "<$gfffile") || die("could not open input file $gfffile");
my ($newfile,$path,$ext) = fileparse($gfffile,qw(\.gff \.gff3));
my $transcripts_out = new Bio::SeqIO(-format => 'Fasta', '-file' => ">$newfile.transcripts.fa");
$newfile =~ s/\.gff//g;
my $db = Bio::DB::Fasta->new($fa);

my @gff;
my %exon_recs;
my $i;
while (<IN>) {
  chomp;
  next if /^\#/o;
  next if /^$/o;
  my @gff = split("\t",$_);
  next if $gff[2] ne 'exon';
  $gff[8]=~/Parent=([^;]+)/;
  my $id = $1;
  $i = 0 if (!exists $exon_recs{$id});
  $exon_recs{$id}->{$i} = $_;
  $i++;
}

foreach my $transcript (sort keys %exon_recs){
  my @sortedexons;
 # print "$exon_recs{$transcript}->{0}\n";
  my @line = split /\t/, $exon_recs{$transcript}->{0};
  if ($line[6] eq '+'){
    #print "$exon_recs{$transcript}->{0}\n";
    @sortedexons = (sort {$a <=> $b} keys %{$exon_recs{$transcript}});
  #  print "@sortedexons\n";
  }
  else {
    @sortedexons = (sort {$b <=> $a} keys %{$exon_recs{$transcript}});
   # print "@sortedexons\n";
  }
  my $sid = $line[0];
  my $ntseq = '';
  foreach my $exons (@sortedexons){
    #print "$exon_recs{$transcript}->{$exons}\n";
    my @line = split /\t/, $exon_recs{$transcript}->{$exons};
    #print "$sid\t$line[3]\t$line[4]\t$line[6]\n";
    $ntseq .= SeqOp::get_seq_BioDBFasta($db, $sid, $line[3], $line[4], $line[6]);    
  } 
  my $exonobj  = Bio::Seq->new( 
                              -seq => $ntseq,
                              -id  => $transcript,
                              -alphabet   => 'dna',
                           );
  $transcripts_out->write_seq($exonobj);

}

