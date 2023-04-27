#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename qw( fileparse );

my $proj = 'EVM';
my $annotation_version = '1A';
my $src = 'CNAG';


# Need to add GFF2 FITS (geneid-like) output
my $version = 3;
my $source = 0;
my $removestop = 0;
my $addstop = 0;
my $txtag = 'transcript';
my $mrna = 0;
my $ov = 3;
my $fits = 0;
my $fl = 0;
my $out = $$;
my $fname = '';
GetOptions(
	   'f|i|a:s'        => \$fname,
	   'v:s' => \$annotation_version,
	   'base|project|species|prefix:s' => \$proj
	  );
my ($base,$path,$ext) = fileparse($fname,qw(\.gff3 \.gff));
die "$fname not a gff3 file" if !$ext;
my %transcripts;
my %tx_gene;
my %genes;
my %tx;
my $cluster = 0;
my $TRANS = '';
my %gene_cluster;
my %cluster_genes;
my %starts;
my %ends;
my %cluster_start = undef;
my %cluster_end = undef;
my %cluster_seq = undef;
my %clustergff;
my $check_start_stop_codon = 0;
#file is well sorted with gene records before transcripts before CDS and exon
#protein sequences follow 
open IN, "<$fname" or die "Couldn't open $fname: $!\n";
my $original = 1;
my $count = 0;
print STDERR "Slurping up GFF records\n";
while (<IN>) {
  chomp;
  if ( m/^[# ]/ || $_!~/\S/) {
    next;
  }
  my @gtf = split "\t",$_;
  next if $gtf[2] !~/CDS|exon|gene|mRNA|transcript/i;  
  if ($gtf[2] eq "gene") {
    #ignore gene records for now
  } elsif ($gtf[2] =~/(mRNA|transcript)/) {
    ### start new cluster for every transcript
    if ($gtf[8]=~m/ID=([^;]+)/) {
      #print STDERR "($1)\n";
      $genes{$1}=\@gtf;
      $cluster++;
      $TRANS=$1;
      push(@{$cluster_genes{$cluster}},$TRANS);
      $cluster_seq{$cluster}=$gtf[0];
      $cluster_start{$cluster}=$gtf[3];
      $cluster_end{$cluster}=$gtf[4];
      $clustergff{$cluster}=\@gtf;
      print STDERR "." if !($count++ % 1000);
    }
    $gtf[2] =$txtag;
    my $gene =0;
    my $id =0;
    if ($gtf[8]=~m/Parent=([^;]+)/) {
      $gene = $1;
    }
    if ($gtf[8]=~m/ID=([^;]+)/) {
      $id = $1;
    }
    if ($id) {
      $tx{$id}->{gff}=\@gtf;
      $tx{$id}->{id}=$id;
      $tx{$id}->{gene}=$gene;
      $tx{$id}->{original}=$original;
      $gene_cluster{$id}=$cluster;
      $tx{$id}->{cluster}=$cluster;
    }
  } elsif ($gtf[2] =~/(exon|CDS)/i) { 
    if ($gtf[8]=~m/Parent=([^;]+)/) {
      my $txid=$1;
      push @{$tx{$txid}->{exons}},\@gtf if $gtf[2] eq "exon";
      push @{$tx{$txid}->{cds}},\@gtf if $gtf[2] eq "CDS";
      if ($gtf[8]=~m/partial[^=]*=true/) {
	$tx{$txid}->{partial}=1;
      } 
      if (exists($starts{$gtf[0]."-".$gtf[3].$gtf[6]})) {
	if ($cluster ne $gene_cluster{$starts{$gtf[0]."-".$gtf[3].$gtf[6]}}) {
	  #print STDERR $starts{$gtf[0]."-".$gtf[3]}," with start: ",$gtf[0]."-".$gtf[3]," was cluster ",$gene_cluster{$starts{$gtf[0]."-".$gtf[3]}}," and is now cluster $cluster\n";
	}
	my $previous_cluster = $gene_cluster{$starts{$gtf[0]."-".$gtf[3].$gtf[6]}};
       	if (exists $cluster_genes{$previous_cluster}) {
	  foreach my $g (@{$cluster_genes{$previous_cluster}}) {
	    #print STDERR "$g\t$previous_cluster\t$cluster\t$genes{$g}->[0]\t$gtf[0]\t$gtf[2]\t$gtf[3]\t$gtf[4]\t$gtf[8]\n";
	    #die "$cluster_seq{$cluster} ne $gtf[0]\n" if $cluster_seq{$cluster} ne $gtf[0];
	    if ($gene_cluster{$g} eq $previous_cluster) {
	      $gene_cluster{$g}=$cluster;
	      my $gene_exists=0;
	      foreach my $cg (@{$cluster_genes{$cluster}}) {
		if ($g eq $cg) {
		  $gene_exists=1;last;
		}
	      }
	      push @{$cluster_genes{$cluster}},$g if !$gene_exists;
	      $cluster_start{$cluster}=$genes{$g}->[3] if $genes{$g}->[3] < $cluster_start{$cluster};
	      $cluster_end{$cluster}=$genes{$g}->[4] if $genes{$g}->[4] > $cluster_end{$cluster};
	    }
	  }
	}
      } else {
	$starts{$gtf[0]."-".$gtf[3].$gtf[6]}=$txid; #$tx{$txid}->{gene};
      }
      if (exists($ends{$gtf[0]."-".$gtf[4].$gtf[6]})) {
	my $previous_cluster = $gene_cluster{$ends{$gtf[0]."-".$gtf[4].$gtf[6]}};
       	if (exists $cluster_genes{$previous_cluster}) {
	  foreach my $g (@{$cluster_genes{$previous_cluster}}) {
	    if ($gene_cluster{$g} eq $previous_cluster) {
	      $gene_cluster{$g}=$cluster;
	      my $gene_exists=0;
	      foreach my $cg (@{$cluster_genes{$cluster}}) {
		if ($g eq $cg) {
		  $gene_exists=1;last;
		}
	      }
	      push @{$cluster_genes{$cluster}},$g if !$gene_exists;
	      $cluster_start{$cluster}=$genes{$g}->[3] if $genes{$g}->[3] < $cluster_start{$cluster};
	      $cluster_end{$cluster}=$genes{$g}->[4] if $genes{$g}->[4] > $cluster_end{$cluster};
	    }
	  }
	}
      } else {
	$ends{$gtf[0]."-".$gtf[4].$gtf[6]}=$txid; #$tx{$txid}->{gene};
      }
    }
  }  
}
close IN;
print STDERR "\n$cluster clusters\n";
print STDERR "\nSorting transcripts...\n";
my @unsorted;
foreach my $t (keys %tx) {
  #my @CDS = sort {$a->start <=> $b->start} @{$t->{CDS}};
  #$t->{CDS}=\@CDS;
  push @unsorted, $tx{$t};
}
my @sorted_transcripts = sort tsort @unsorted;

my %seen_genes;
my %seen_cds;
my %seen_structure;
my $gene_count = 1;
my %tcount = 1;
my %ccount = 1;
my %pcount = 1;
my $newgenename = '';
open GFF,">$$.$proj.$annotation_version.gff3" or die "Couldn't open $proj.$annotation_version.gff3: $!";
print GFF "##gff-version 3\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf GFF "# date: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
my $first = 1;
foreach my $transcript (@sorted_transcripts) {
  if (!(exists($transcript->{cds}))) {
    print STDERR $transcript->{id}," has no CDS\n"; next;
  }
  my $canon = 0;
  my @tgff = @{$transcript->{gff}};
  my $strand = $tgff[6];
  my $tid = $transcript->{id};
  my $gene_id = $gene_cluster{$tid};
  $newgenename = sprintf("%s%s%06d",$proj,$annotation_version,$gene_id);
  if (!exists $seen_genes{$gene_id}) {
    print GFF "###\n" if ($ov==3 && !$first);
    $first=0;
    print GFF (join("\t",($tgff[0],$src,'gene',$cluster_start{$gene_id},$cluster_end{$gene_id},'.',$tgff[6],'.',"ID=$newgenename")),"\n");
    $seen_genes{$gene_id}++;
    $tcount{$newgenename} = 1;
    $ccount{$newgenename} = 0;
    $pcount{$newgenename} = 0;
  }
  my $uniq_str = '';
  foreach my $gr (sort gffsort @{$transcript->{cds}}) {
    $uniq_str .=join(':',('c',$gr->[0],$gr->[3],$gr->[4],$gr->[6]));
  }
  $uniq_str .=":";
  foreach my $gr (sort gffsort @{$transcript->{exons}}) {
    $uniq_str .=join(':',('e',$gr->[0],$gr->[3],$gr->[4],$gr->[6]));
  }
  ### transcript sorted by PASA_UPDATE's followed by ORIGINAL
  if (! exists $seen_structure{$uniq_str}) {
    #$seen_structure{$uniq_str} = $transcript;
    $seen_structure{$uniq_str}++;
  } else {
    print STDERR "$uniq_str seen already -- skipping\n";
    next;
  }
  my $newcdsid = $newgenename.'C'.$tcount{$newgenename};
  my $newtid = $newgenename.'T'.$tcount{$newgenename}++;
  $transcript->{gff}->[1]=$src;
  $transcript->{gff}->[8]=~s/ID=[^;]+/ID=$newtid/;
  $transcript->{gff}->[8]=~s/Parent=[^;]+/Parent=$newgenename/;
  $transcript->{gff}->[8]=~s/;Name=[^;]+/;Name=$newtid/;
  #print GFF (join("\t",@{$transcript->{gff}}),"\n");
 
  my @cds;
  my $count = 1;
  my $ecount = 1;
  my $numexon = 0;
  my $numcds = 0;
  my $string = '';
  my $cdslen = 0;
  foreach my $gffrec (sort gffsort @{$transcript->{cds}}) {
    $cdslen+=$gffrec->[4]-$gffrec->[3]+1;
    $numcds++;
    $string .=join(':',($gffrec->[0],$gffrec->[3],$gffrec->[6],$gffrec->[4],$gffrec->[6]));
  }
  my $newprotid = $newgenename.'P'.$pcount{$newgenename};
  if (!exists $seen_cds{$string}) {
    $pcount{$newgenename}++;
    $newprotid = $newgenename.'P'.$pcount{$newgenename};
    $seen_cds{$string}=$newprotid;
  } else {
    $newprotid = $seen_cds{$string};
  }
  $transcript->{gff}->[8]=~s/;?product=[^;]+//;
  $transcript->{gff}->[8]=~s/;$//;
  $transcript->{gff}->[8].=";product=$newprotid";
  if (exists($transcript->{partial}) && $transcript->{partial}) {
    $transcript->{gff}->[8].=";partial_cds=true";
  }
  print GFF (join("\t",@{$transcript->{gff}}),"\n");
  if ($strand eq '+') {
    foreach my $gffrec (sort gffsort @{$transcript->{exons}}) {
      my $newexonid = $newtid.'.exon'.$ecount;
      $gffrec->[8]=~s/ID=[^;]+/ID=$newexonid/;
      $gffrec->[8]=~s/Parent=[^;]+/Parent=$newtid/;
      $gffrec->[1]=$src;
      $gffrec->[8]=~s/;?Name=[^;]+//;
      $gffrec->[8].=";Name=$newtid";
      print GFF (join("\t",@{$gffrec}),"\n");
      $ecount++;
    }
  } else {
    my $ecount = scalar @{$transcript->{exons}};
    foreach my $gffrec (sort gffsort @{$transcript->{exons}}) {
      my $newexonid = $newtid.'.exon'.$ecount;
      $gffrec->[8]=~s/ID=[^;]+/ID=$newexonid/;
      $gffrec->[8]=~s/Parent=[^;]+/Parent=$newtid/;
      $gffrec->[1]=$src;
      $gffrec->[8]=~s/;?Name=[^;]+//;
      $gffrec->[8].=";Name=$newtid";
      print GFF (join("\t",@{$gffrec}),"\n");
      $ecount--;
    }
  }
  #print all CDS records
  $canon = 0; 
  my $pstart = 1;
  my $sum = 0;
  my $first = 1;
  my $numcds = scalar  @{$transcript->{cds}};
  my $cdscount = 0;
  my @sortedcds = sort gffsort @{$transcript->{cds}};
    
  # need to match up protein length with cds length. adjust first exon length by phase or remainder then keep summing exon lengths.
  my $phase = $sortedcds[0]->[7];
  my $rmd = ($cdslen - $phase)%3;
  $first = 1;
  foreach my $gffrec (sort gffsort @{$transcript->{cds}}) {
    $gffrec->[8]=~s/partial_gene=true/partial_cds=true/;
    $gffrec->[8]=~s/ID=[^;]+/ID=$newcdsid/;
    $gffrec->[8]=~s/Parent=[^;]+/Parent=$newtid/;
    $gffrec->[8]=~s/;Target=[^;]+//;
    $gffrec->[8]=~s/;?Name=[^;]+//;
    $gffrec->[8].=";Name=$newcdsid";
    #    $gffrec->[8].=";Target=$newprotid $pstart $pend";
    $gffrec->[1]=$src;
    print GFF (join("\t",@{$gffrec}),"\n");
  }
  $ccount{$newgenename}++;
}
close GFF;
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
