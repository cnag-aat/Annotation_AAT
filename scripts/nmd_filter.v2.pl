#!/usr/bin/env perl
use strict;
#use Getopt::Long;
#use File::Basename qw( fileparse );

my $nmd_dist=50;
my $nmd_dist_soft = 150;
my %transcript=();
my %genes;
my %geneprinted=();
my $tid='';
print STDERR "Tagging likely NMD targets. Transcripts must be grouped by gene. Gene record first followed by transcript, exons, cds and followed by ###, even last line.\n";
while (<>) {
  chomp;
  if ( m/^###/) {
    #next if !defined $transcript{$tid}->{gff};
    foreach my $t (sort keys %transcript) {
      print STDERR "Non-coding transcript?\n" if ! exists $transcript{$t}->{cds};
      my @sortedcds = sort gffsort @{$transcript{$t}->{cds}};
      my $nmd=0;
      my $nmdsoft=0;
      my $dist=0;
      if ($transcript{$t}->{strand} eq "+") {
	my $stop=$sortedcds[-1]->[4];
	foreach my $e (sort gffsort @{$transcript{$t}->{exons}}) {
	  $nmd=1 if $dist > $nmd_dist;
	  $nmdsoft=1 if $dist > $nmd_dist_soft;
	  if ($e->[4]>$stop) {
	    if ($e->[3]>$stop) {
	      $dist+=($e->[4] - $e->[3] + 1);
	    } else {
	      $dist+=($e->[4] - $stop);
	    }
	  }
	}
      } else {			# - strand
	my $stop=$sortedcds[0]->[3];
	foreach my $e (sort reverse_gffsort @{$transcript{$t}->{exons}}) {
	  $nmd=1 if $dist > $nmd_dist; 
	  $nmdsoft=1 if $dist > $nmd_dist_soft;
	  if ($e->[3]<$stop) {
	    if ($e->[4]<$stop) {
	      $dist+=($e->[4] - $e->[3] + 1);
	    } else {
	      $dist+=($stop - $e->[3]);
	    }
	  }
	}
      }
      if (! exists $geneprinted{$transcript{$t}->{gid}}) {
	print $transcript{$t}->{gene};
	$geneprinted{$transcript{$t}->{gid}}++;
      }
      if($nmd || $nmdsoft){
	$transcript{$t}->{tgff}=~ s/transcript_biotype=protein_coding/transcript_biotype=non_sense_mediated_decay/;
        $transcript{$t}->{gff} =~ s/transcript_biotype=protein_coding/transcript_biotype=non_sense_mediated_decay/g;
      }
      print  $transcript{$t}->{tgff};
      print  $transcript{$t}->{gff};
      delete $transcript{$t};
    }
    print "###\n";
  } elsif ( m/^#/) {
    print "$_\n";
  } else {
    my @gff = split "\t",$_;
    next if $gff[2] !~/CDS|exon|gene|mRNA|transcript/i;
    
    #print STDERR $gff[8],"\n";
    if ($gff[2] eq "gene") {
      chomp $gff[8];
      if ($gff[8]=~m/ID=([^;]+)/) {
	$genes{$1}=$_ . ";gene_type=protein_coding\n";
	#print STDERR "$1\n$_";
      }
    } elsif ($gff[2] =~/(mRNA|transcript)/) {
      my $gene;
      chomp $gff[8];
      if ($gff[8]=~m/Parent=([^;]+)/) {
	$gene = $1;
	die "Gene record must come before transcripts for that gene\n" if ! exists $genes{$gene};
      }
      if ($gff[8]=~m/ID=([^;]+)/) {
	$tid = $1;
      }
      if ($tid && $gene) {
	$transcript{$tid}->{gid}=$gene;
	$transcript{$tid}->{gene}=$genes{$gene};
	$transcript{$tid}->{tgff}=$_ . ";gene_type=protein_coding;transcript_biotype=protein_coding\n";
	$transcript{$tid}->{strand}=$gff[6];
      }
    } elsif ($gff[2] =~/CDS/i) {
    
      if ($gff[8]=~m/Parent=([^;]+)/) {
	push @{$transcript{$1}->{cds}},\@gff;
	$transcript{$tid}->{gff}.=$_ . ";gene_type=protein_coding;transcript_biotype=protein_coding\n";
      }
    } elsif ($gff[2] =~/exon/i) {
      if ($gff[8]=~m/Parent=([^;]+)/) {
	push @{$transcript{$1}->{exons}},\@gff;
	$transcript{$tid}->{gff}.=$_ . ";gene_type=protein_coding;transcript_biotype=protein_coding\n";
      }
    }
  }
}


print STDERR "done!\n";
sub tsort
  {
    $a->{gff}->[0] cmp $b->{gff}->[0]
      ||
	$a->{gff}->[3] <=> $b->{gff}->[3]
	  ||
	    $a->{gff}->[4] <=> $b->{gff}->[4]
	      ||
		$a->{gff}->[6] <=> $b->{gff}->[6]
		  ||
		    $a->{original} <=> $b->{original}
		  
	      }
sub gffsort
  {
    $a->[0] cmp $b->[0]
      ||
	$a->[3] <=> $b->[3]
	  ||
	    $a->[4] <=> $b->[4]
	      ||
		$a->[6] <=> $b->[6]
	      }

sub reverse_gffsort
  {
    $b->[0] cmp $a->[0]
      ||
	$b->[3] <=> $a->[3]
	  ||
	    $b->[4] <=> $a->[4]
	      ||
		$b->[6] <=> $a->[6]
	      }
