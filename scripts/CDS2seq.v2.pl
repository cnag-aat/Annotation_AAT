#!/usr/bin/env perl

############################################################
#
# script to extract CDS nt and protein sequence from CDS records in GFF3 format
# 
# Author: T. Alioto
# Date: 20130905
############################################################

use strict;
use lib "/home/groups/assembly/talioto/myperlmods/";
#use lib "/apps/BIOPERL/1.6.1/lib/perl5";
#use lib "/scratch/project/devel/aateam/perlmods";
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Basename qw( fileparse );
use Getopt::Long;
use SeqOp;

my $gfffile =0;	
my $fa= 0;

my $code = 1;
GetOptions(
	   'i|g|file|gff:s'		=> \$gfffile,
	   'f|s|seq:s'		=> \$fa,
	   'c|code:s'		=> \$code
	  );

   
############################################################
#$fragseq = SeqOp::get_seq_BioDBFasta($db, $sid,$forstart,$revend,'+');
open (IN, "<$gfffile") || die("could not open input file $gfffile");
my ($newfile,$path,$ext) = fileparse($gfffile,qw(\.gff \.gff3));
my $nt_out = new Bio::SeqIO(-format => 'Fasta', '-file' => ">$newfile.cds.fa");
my $aa_out = new Bio::SeqIO(-format => 'Fasta', '-file' => ">$newfile.pep.fa");
$newfile =~ s/\.gff//g;
#open (CDS, ">$newfile".".CDS.fa");
#if ($pretty){ 
#	open PRETTY,">$newfile"."_splice_jxns.txt";
#	if ($header) {
#		print PRETTY "$header\n";
#	}
#}

#if ($print_intron_seq) {
#	open (ISEQ, ">$newfile"."_full_intron_seq.fa");
#}
my $db = Bio::DB::Fasta->new($fa);
my $seq = $db->get_Seq_by_id('Pcu23_ss01');
#print "$seq\n";
my @gff;
my %cds_recs;
my %p;
while (my $l = <IN>) {
  if ($l =~ m/^[^# ]/) { #Skip gff comments and metainfo and blank lines
    chomp $l;		
    my @gff = split("\t",$l);
    if ($gff[2]=~/transcript|mRNA/){
      $gff[8]=~m/product=([^;]+)/;
      my $protein = $1;
      $gff[8]=~m/ID=([^;]+)/;
      my $tid=$1;
      $p{$tid}=$protein;
    }
    next if $gff[2] ne 'CDS';
    $gff[8]=~/Parent=([^;]+)/;
    my $id = $1;
    #if ( $gff[8]=~/Target=([^; ]+)/){
      #print STDERR "$1\n";
      #$id .= "|$1";
    #}
    push @{$cds_recs{$id}},\@gff;
  }
}
my %seen = ();
foreach my $cdsr (sort keys %cds_recs){
  print "$cdsr\n";
  if(exists $seen{$p{$cdsr}}){
    next;
  }else{
    $seen{$p{$cdsr}}++;
  }
  my @sortedexons;
  if ($cds_recs{$cdsr}->[0]->[6] eq '+'){
    @sortedexons = sort {$a->[3]<=>$b->[3]} @{$cds_recs{$cdsr}};
  }else{
    @sortedexons = sort {$b->[3]<=>$a->[3]} @{$cds_recs{$cdsr}};
  }
  my $sid=$sortedexons[0]->[0];
  my $ntseq = '';
  foreach my $cds (@sortedexons){
    $ntseq .= SeqOp::get_seq_BioDBFasta($db, $sid,$cds->[3],$cds->[4],$cds->[6]);
  }
  #print "$sid\t$ntseq\n";
  my $description = '';
  if ($ntseq !~ /^ATG.*(TAA|TGA|TAG)$/ || $sortedexons[0]->[7]){$description = "partial_cds";}
  my $cdsobj  = Bio::Seq->new( 
			      -seq => $ntseq,
			      -description => $description,
			    -id  => $p{$cdsr},
			    -alphabet 	=> 'dna',
			   );
  $nt_out->write_seq($cdsobj);
  my $protcdsobj = $cdsobj->translate(undef,undef,$sortedexons[0]->[7],$code);
  $aa_out->write_seq($protcdsobj);
}

############################################################
