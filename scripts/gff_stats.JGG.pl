#!/usr/bin/env perl
use strict;
#use lib "/home/devel/talioto/myperlmods";
use lib "/home/groups/assembly/talioto/myperlmods/";
use Getopt::Long;
use SeqOp;
use Bio::DB::Fasta;
use File::Basename;
my $genome_fasta = '';
my $annotation = '';
 GetOptions(
	    'f|s|g:s'=>\$genome_fasta,
	    'a|i:s'=>\$annotation
 	   );
my $SCRATCH = $ENV{'TMPDIR'};
#`mkdir -p $SCRATCH/$uid`;
my ($abase,$apath,$aext) = fileparse($annotation,qw(\.gff \.gff3));
my ($base,$path,$ext) = fileparse($genome_fasta,qw(\.fa \.fasta \.mfa));
if (!-e "$SCRATCH/$base$ext") {
  `cp $genome_fasta $SCRATCH`;
}
my $genome_mfa_orig = $genome_fasta;
$genome_fasta = "$SCRATCH/$base$ext";

#print STDERR "$genome_fasta\n";
if (!-e "$genome_fasta.index") {
  print STDERR "indexGenome.pl -f $genome_fasta\n";
  `indexGenome.pl -f $genome_fasta`;
}
die "$genome_fasta.index does not exist!\n" if !-e "$genome_fasta.index";
my $db = Bio::DB::Fasta->new($genome_fasta);
my @seqids = sort $db->ids;
die if ! @seqids;

my $genome_size = 0;
foreach my $id (@seqids) {
  $genome_size+=$db->length($id);
}


my %transcripts;
my %genes;
my @exons;
my %uniqexons;
my @cds;
my %uniqcexons;
my %pepreg;
my %uniqintrons;
my $exonseq = '';
my $cdsseq = '';
my $intronseq = '';
my @totUTRlen;
my @UTR5;
my @UTR3;
my $multiexonic = 0;
my $numt = 0;
my $numex=0;
my $numint = 0;
open ANN,"<$annotation" or die "Couldn't open $annotation: $!\n";
print STDERR "Reading $annotation\n";
while(<ANN>){
    next if m/^[# ]/;
    chomp;
    my @gff = split "\t",$_;
    next if $gff[2] !~/CDS|First|Internal|Terminal|initial|exon|gene|transcript|mRNA|start_codon|stop_codon|Intron/i;
    my %att;
	#print STDERR $gtf[8],"\n";
    if ($gff[2] eq "gene"){
      if ($gff[8]=~m/ID=([^;]+)/){
	$genes{$1}->{gff}=\@gff;
	$genes{$1}->{length}=$gff[4]-$gff[3]+1;
      }else{
	die "gene has no ID\n";
      }
    }elsif($gff[2] =~ /mRNA|transcript/){
      my $t = '';
      my $g = '';
      $numt++;
      if ($gff[8]=~m/ID=([^;]+)/){
	$t=$1;
	$transcripts{$t}->{gff}=\@gff;
	$transcripts{$t}->{strand}=$gff[6];
        $transcripts{$t}->{length}=0;
        $pepreg{$t}->{length}=0;
      }else{
	die "transcript has no ID\n";
      }
    }elsif ($gff[2] =~/exon/i) {
      $numex++;
      $uniqexons{$gff[0],$gff[3],$gff[4]}->{gff}=\@gff;
      $uniqexons{$gff[0],$gff[3],$gff[4]}->{length}=$gff[4]-$gff[3]+1;
      my $t = '';
      if ($gff[8]=~m/Parent=([^;]+)/){
	$t=$1;
        $transcripts{$t}->{length}+=$gff[4]-$gff[3]+1;
        push @{$transcripts{$t}->{exons}},{gff=>\@gff};
      }else{
	die "exon has no Parent\n";
      }
    }
    elsif ($gff[2] =~/CDS|First|Internal|Terminal|Single/i) {
      $uniqcexons{$gff[0],$gff[3],$gff[4]}->{gff}=\@gff;
      $uniqcexons{$gff[0],$gff[3],$gff[4]}->{length}=$gff[4]-$gff[3]+1;
      my $t = '';
      if ($gff[8]=~m/Parent=([^;]+)/){
	$t=$1;
	#my $g = $transcripts{$t}->{gene};
        $pepreg{$t}->{length}+=$gff[4]-$gff[3]+1;
        push @{$transcripts{$t}->{cds}},{gff=>\@gff};
      }else{
	die "exon has no Parent\n";
      }
    }
}
 #my @unsorted;
print STDERR "Parsing transcripts\n";
foreach my $t (keys %transcripts){
  my @e = sort {$a->{gff}->[3] <=> $b->{gff}->[3]} @{$transcripts{$t}->{exons}};
  my $first_ex = 1;
  my $pEND =0;
  $multiexonic++ if scalar @e > 1;
  foreach my $ex (@e){
    my $seq = $ex->{gff}->[0]; 
    my $start = $ex->{gff}->[3];
    my $end = $ex->{gff}->[4];
    if ($first_ex){
      $first_ex=0;
    }else{
      $uniqintrons{$seq,$start,$end}->{length} = $start-$pEND+1;
      $intronseq  .= SeqOp::get_seq_BioDBFasta($db,  $ex->{gff}->[0],$pEND,$start,'+');
      $numint++;
    }
    $pEND=$ex->{gff}->[4];
  }
  my @coding_exons = sort {$a->{gff}->[3] <=> $b->{gff}->[3]} @{$transcripts{$t}->{cds}};
  my $UTR5len = 0;
  my $UTR3len = 0;
  if ($transcripts{$t}->{strand} eq '+'){
    my $start = $coding_exons[0]->{gff}->[3];
    my $stop = $coding_exons[-1]->{gff}->[4];
    foreach my $ex (@e){
      my $seq = $ex->{gff}->[0];
      my $estart =$ex->{gff}->[3];
      my $eend = $ex->{gff}->[4];
      if ($eend<$start){
        $UTR5len += $eend - $estart + 1;
      }elsif($estart<$start){
        $UTR5len += $start - $estart;
      }elsif($estart>$stop){
        $UTR3len += $eend - $estart + 1;
      }elsif($eend>$stop){
        $UTR3len += $eend - $stop;
      }
    }
  }else{
     my $stop = $coding_exons[0]->{gff}->[3];
     my $start = $coding_exons[-1]->{gff}->[4];
     foreach my $ex (@e){
       my $seq = $ex->{gff}->[0];
       my $estart =$ex->{gff}->[3];
       my $eend = $ex->{gff}->[4];
       if ($eend<$stop){
        $UTR3len += $eend - $estart + 1;
       }elsif($estart<$stop){
        $UTR3len += $stop - $estart;
       }elsif($estart>$start){
        $UTR5len += $eend - $estart + 1;
       }elsif($eend>$start){
        $UTR5len += $eend - $start;
       }
     }
  }
  push @totUTRlen, {length=>$UTR5len + $UTR3len} if ($UTR5len + $UTR3len);
  push @UTR5,{length=>$UTR5len} if $UTR5len;
  push @UTR3,{length=>$UTR3len} if $UTR3len;
}

close ANN;
############ BEGIN STATS #############
#### GENOME ####
print STDERR "Computing statistics\n";
print "Genome size (bp)\t$genome_size\n";

#### GENES ####
my @genes;
foreach my $key (sort keys %genes){
  #print STDERR "$key\n";
  push @genes,$genes{$key};
}
my %genestats = num_mean_median(\@genes);
print "Number of genes\t$genestats{n}\n";
print "Mean gene length\t$genestats{mean}\n";
print "Median gene length\t$genestats{median}\n";
print "Total genic length\t$genestats{total}\n";
print "Gene length range\t$genestats{min}..$genestats{max}\n";
print "Gene density (gene/Mb)\t",$genestats{n}/($genome_size/1000000),"\n";

#### EXONS ####
my @exons; 
foreach my $key (sort keys %uniqexons){
  #print STDERR "$key\n";
  push @exons,$uniqexons{$key};
}
my %exonstats = num_mean_median(\@exons);
my $nexons = keys %uniqexons;
print "Number of exons\t$nexons\n";
print "Mean exon length\t$exonstats{mean}\n";
print "Median exon length\t$exonstats{median}\n";
print "Total exon length\t$exonstats{total}\n";
print "Exon length range\t$exonstats{min}..$exonstats{max}\n";
#### Coding exons ####
my @coding_exons;
foreach my $key (sort keys %uniqcexons){
  #print STDERR "$key\n";
  push @coding_exons,$uniqcexons{$key};
  $cdsseq  .= SeqOp::get_seq_BioDBFasta($db, $uniqcexons{$key}->{gff}->[0],$uniqcexons{$key}->{gff}->[3],$uniqcexons{$key}->{gff}->[4],'+');
}
#print "$cdsseq\n";
my %cds_seg_stats = num_mean_median(\@coding_exons);
print "Number of coding exons\t$cds_seg_stats{n}\n";
print "Mean coding exon length\t$cds_seg_stats{mean}\n";
print "Median coding exon length\t$cds_seg_stats{median}\n";
print "Total coding exon length\t$cds_seg_stats{total}\n";
print "Coding exon length range\t$cds_seg_stats{min}..$cds_seg_stats{max}\n";
#### INTRONS ####
my @introns;
foreach my $key (sort keys %uniqintrons){
  #print STDERR "$key\n";
  push @introns,$uniqintrons{$key};
}
my %intronstats = num_mean_median(\@introns);
print "Number of introns\t$intronstats{n}\n";
print "Mean intron length\t$intronstats{mean}\n";
print "Median intron length\t$intronstats{median}\n";
print "Intron length range\t$intronstats{min}..$intronstats{max}\n";
#my $igc = gcpct($intronseq);
#chomp $igc;
#print "Intron GC content\t$igc\n";
#### Transcripts ####
my @txlengths; 
foreach my $key (sort keys %transcripts){
  #print STDERR "$key\n";
  push @txlengths,$transcripts{$key};
}
my %tstats = num_mean_median(\@txlengths);
print "Number of transcripts\t$tstats{n}\n";
print "Mean transcript length\t$tstats{mean}\n";
print "Median transcript length\t$tstats{median}\n";
print "Transcript length range\t$tstats{min}..$tstats{max}\n";
#### CDS ####
my @cdslengths; 
foreach my $key (sort keys %pepreg){
  #print STDERR "$key\n";
  push @cdslengths,$pepreg{$key};
}
my %cdslstats = num_mean_median(\@cdslengths);
print "Number of cds\t$cdslstats{n}\n";
print "Mean cds length\t$cdslstats{mean}\n";
print "Median cds length\t$cdslstats{median}\n";
print "Total cds length\t$cdslstats{total}\n";
print "CDS length range\t$cdslstats{min}..$cdslstats{max}\n";
#my $cgc = gcpct($cdsseq);
#chomp $cgc;
#print "Coding GC content\t$cgc\n";

#### UTRs ####
my %utrstats = num_mean_median(\@totUTRlen);
print "Number of transcripts with UTRs\t$utrstats{n}\n";
print "Mean 5' + 3' UTR length\t$utrstats{mean}\n";
print "Median 5' + 3' UTR length\t$utrstats{median}\n";
print "Total 5' + 3' UTR length\t$utrstats{total}\n";
print "5' + 3' UTR range\t$utrstats{min}..$utrstats{max}\n";
my %utr5stats = num_mean_median(\@UTR5);
print "Number of transcripts with 5' UTRs\t$utr5stats{n}\n";
print "Mean 5' UTR length\t$utr5stats{mean}\n";
print "Median 5' UTR length\t$utr5stats{median}\n";
print "Total 5' UTR length\t$utr5stats{total}\n";
print "5' UTR range\t$utr5stats{min}..$utr5stats{max}\n";
my %utr3stats = num_mean_median(\@UTR3);
print "Number of transcripts with 3' UTRs\t$utr3stats{n}\n";
print "Mean 3' UTR length\t$utr3stats{mean}\n";
print "Median 3' UTR length\t$utr3stats{median}\n";
print "Total 3' UTR length\t$utr3stats{total}\n";
print "3' UTR range\t$utr3stats{min}..$utr3stats{max}\n";

print "\n";
print "Exons per transcript\t",$numex/$numt,"\n";
print "Introns per transcript\t",$numint/$numt,"\n";
print "Transcripts per gene\t",$numt/$genestats{n},"\n";
print "Multi-exonic transcripts\t",$multiexonic/$numt,"\n";

sub num_mean_median{
  my $array = shift;
  my $tot=0;
  my $n = scalar @$array;
  return undef if $n < 1;
  my @sorted = sort {$a->{length} <=> $b->{length}} @$array;
  my $median = $sorted[int($n/2)]->{length};
  my $min = $sorted[0]->{length};
  my $max = $sorted[-1]->{length};
  foreach my $hashref (@sorted){
    $tot+=$hashref->{length};
  }
  my $mean = $tot/$n;
  return (n=>$n,mean=>$mean,median=>$median,total=>$tot,min=>$min,max=>$max);
}
system("exon_intron_lengthdist.R $abase");
print STDERR "done!\n";
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

sub gcpct{
  my %gc; my $t=0;
  my $seq = shift;
  chomp $seq;
  while ($seq=~/([ACGTacgt])/g){
    $gc{uc($1)}++;
    $t++;
  }


  foreach my $k (sort keys %gc){
    #printf "%s\t%2.3f\n",$k,$gc{$k}/$t;
  }
  return sprintf("%2.2f%%\n",($gc{G}/$t + $gc{C}/$t)*100);
  #return $t;
}
