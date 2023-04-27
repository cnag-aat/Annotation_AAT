#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($chunkfile, $species, $alternatives, $sample, $gff3, $noInFrameStop, $uniqueGeneId, $maxtracks, $strand, $singlestrand, $intronlength, $extrinsic_file,
$additional);

GetOptions(
           'f:s'             => \$chunkfile,
           's:s'             => \$species,
           'alt:s'           => \$alternatives,
           'sample:s'        => \$sample,
           'gff:s'           => \$gff3,
           'iF:s'            => \$noInFrameStop,
           'uniq:s'          => \$uniqueGeneId,
           'max:s'           => \$maxtracks,
           'str:s'           => \$strand,
           'ss:s'            => \$singlestrand, 
           'int:s'           => \$intronlength,
           'ef:s'            => \$extrinsic_file,
           'ad:s'            => \$additional
          );

my $hints = "hints.sorted.gff.gz";

my %scaff;

while (<STDIN>) {
    chomp;
    my @line = split /\s+/, $_;
    $scaff{$line[0]} = "" if (!exists $scaff{$line[0]});
    $scaff{$line[0]} .= "$_\n";
}

open FASTA, "<", "$chunkfile";
my $id; 
my %keep; 
while (<FASTA>) {
    chomp;
    if ( $_ =~ m/\>/ ){
      my @ids = split '>', $_;
      $id = $ids[1];
    #  close SEQFILE if (SEQFILE);
      open HINTFILE, ">", "$id.hints.gff";
      print HINTFILE "$scaff{$id}";
      close HINTFILE;
      open SEQFILE, ">", "$id.masked.fa";
      $keep{$id}++;
    }
    print SEQFILE "$_\n";  
}
close FASTA;

foreach my $seq (keys %keep){
  system "augustus --species=$species --hintsfile=$seq.hints.gff --alternatives-from-sampling=$alternatives --sample=$sample --gff3=$gff3 --noInFrameStop=$noInFrameStop --uniqueGeneId=$uniqueGeneId --maxtracks=$maxtracks --strand=$strand --singlestrand=$singlestrand --min_intron_len=$intronlength --extrinsicCfgFile=$extrinsic_file $additional $seq.masked.fa > $seq.augustus_introns.gff3";
}
