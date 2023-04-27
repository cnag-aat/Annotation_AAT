#!/usr/bin/env perl
use strict;
use Getopt::Long;

# Need to add GFF2 FITS (geneid-like) output
my $version = 2.5;
my $source = 0;
my $removestop = 0;
my $addstop = 1;
my $txtag = 'transcript';
my $mrna = 0;
my $ov = 3;
my $fits = 0;
my $fl = 0;
GetOptions(
	   'v|version:s'      => \$version,
	   'rmstop|nostop'  => \$removestop,
	   'addstop!'        => \$addstop,
	   'txtag:s'        => \$txtag,
	   'mrna'           => \$mrna,
	   'ov:s'           => \$ov,
	   'fits'           => \$fits,
	   'fl'             => \$fl
	   );

if ($mrna){$txtag = 'mRNA';}
my %transcripts;
my %tx_gene;
my %genes;
my %tx;
my $check_start_stop_codon = 0;
while(<>){
    next if m/^[# ]/;
    chomp;
    my @gtf = split "\t",$_;
    next if $gtf[2] !~/CDS|First|Internal|Terminal|initial|exon|gene|transcript|mRNA|start_codon|stop_codon/i;
    my %att;
    if ($version eq "2.5"){
	#print STDERR $gtf[8],"\n";
	$gtf[8]=~s/;\s*$//;
	my @list = split ";",$gtf[8];
	foreach my $kv (@list) {
	    $kv =~ /(\S+)\s+(\S.*);?/;
	    $att{$1}=$2;
	}
	$gtf[8]='';
	
	if (exists $att{transcript_id} && exists $att{Parent}){delete $att{Parent};}
	foreach my $a (sort keys %att) {
	    $att{$a}=~s/;$//;
	    $att{$a}=~s/^\s*//;
	    $att{$a}=~s/"$//;
	    $att{$a}=~s/^"//;
	    #$gtf[8].="$a=$att{$a};";
	    if ($gtf[2]=~/gene/){
	      if ($a eq 'gene_id') {
		$gtf[8].="ID=$att{$a};";
	      }else{
		$gtf[8].="$a=$att{$a};";
	      }
	      
	    }elsif ($gtf[2]=~/transcript|mRNA/){
	      if ($a eq 'transcript_id') {
		$gtf[8].="ID=$att{$a};";
	      }elsif ($a eq 'gene_id') {
		$gtf[8].="Parent=$att{$a};";
	      }else{
		$gtf[8].="$a=$att{$a};";
	      }
	      
	    }else{
	      if ($a eq 'transcript_id') {
		$gtf[8].="Parent=$att{$a};";
	      }elsif ($a eq 'gene_id') {
		#$gtf[8].="gene=$att{$a};";
	      }else{
		$gtf[8].="$a=$att{$a};";
	      }
	  }
	}
	$gtf[8]=~s/;$//;
	if ($gtf[2] eq "gene" && (exists $att{gene_id} )){
	  $gtf[8]="ID=$att{gene_id}";
	    $genes{$att{gene_id}}=\@gtf;
	}elsif($gtf[2] =~/transcript|mRNA/ && (exists $att{gene_id} )){
	    $tx_gene{$att{transcript_id}}=$att{gene_id};
	}elsif($gtf[2] =~/codon/ && (exists $att{gene_id} )){
	    $check_start_stop_codon = 1;
	}

    }elsif ($version eq "2"){
	#print STDERR $gtf[8],"\n";
	my $grp = $gtf[8];
	$gtf[8]="Parent=$gtf[8];";
	$att{transcript_id}=$grp;
	if ($gtf[2] eq "gene"){
	    $gtf[8]="ID=$grp;";
	    $genes{$grp}=\@gtf;
	}elsif ($gtf[2] =~/First|Internal|Terminal|Single|initial|internal|terminal/i){
	    $gtf[2] = 'CDS';
	}
	
    }elsif($version eq "3"){
	#print STDERR $gtf[8],"\n";
	if ($gtf[2] eq "gene"){
	    if ($gtf[8]=~m/ID=([^;]+)/){
		$genes{$1}=\@gtf;
	    }
	}elsif($gtf[8]=~m/Parent=([^;]+)/){
	    $att{transcript_id}=$1;
	}
    }
    if ($gtf[2] =~/transcript|mRNA/){
	$gtf[2] =$txtag;
	my $gene =0;
	my $id =0;
	if ($gtf[8]=~m/Parent=([^;]+)/){
	    $gene = $1;
	}
	if ($gtf[8]=~m/ID=([^;]+)/){
	    $id = $1;
	}
	if ($id && $gene){
	    $tx_gene{$id}=$gene;
	    #print STDERR $gtf[7],"\n";
	    $tx{$id}=\@gtf;
	    #print STDERR "$id\t$gene\n";
	}
    }elsif ($gtf[2] =~/CDS|First|Internal|Terminal|initial|exon|start_codon|stop_codon/i && exists $att{transcript_id}) {
      if ($gtf[2] =~/CDS|exon/){
	$gtf[8] =~s/ID=([^;]+);?//;
      }
      push @{$transcripts{$att{transcript_id}}},\@gtf;
    }

}
my @unsorted;
foreach my $t (keys %transcripts){
  #my @CDS = sort {$a->start <=> $b->start} @{$t->{CDS}};
  #$t->{CDS}=\@CDS;
  push @unsorted, $transcripts{$t};
}
my @sorted_transcripts = sort gffsort @unsorted;

foreach my $transcript (@sorted_transcripts){
    my @tgff = @{$transcript->[0]};
    $tgff[2]= $txtag;
    $tgff[4]= $transcript->[-1]->[4];
    $tgff[5]='.';
    $tgff[8]=~/ID=([^;]*)/;
    my $originalid = $1;
    #print STDERR "original $originalid\n";
    $tgff[8]=~/Parent=([^;]*)/;
    $tgff[8]="ID=$1";
    $tgff[7]='.';
    if ($ov == 2){$tgff[8]=$1;}
    my $tid = $1;
    #print STDERR "new $tid\n";
   
    if ($ov !=2){
      if (exists $tx_gene{$tid}) {
	#print STDERR "$tid\t".$tx_gene{$tid}."\n";
	if (exists $genes{$tx_gene{$tid}}) {
	  print (join("\t",@{$genes{$tx_gene{$tid}}}),"\n");
	  $tgff[8].=";Parent=$tx_gene{$tid}";
    
	} else {
	  $tgff[8].=";Parent=g_$tid";
	  print (join("\t",($tgff[0],$tgff[1],'gene',$tgff[3],$tgff[4],'.',$tgff[6],'.',"ID=g_$tid")),"\n");
	}
	delete $tx_gene{$tid};
      } else {
	 $tgff[8].=";Parent=g_$tid";
	 print (join("\t",($tgff[0],$tgff[1],'gene',$tgff[3],$tgff[4],'.',$tgff[6],'.',"ID=g_$tid")),"\n");
      }
    }
    my @cds;
    my $count = 1;
    my $numexon = 0;
    my $numcds = 0;
    
    ### Make sure CDS starts and ends with stop. If not it's a partial annotation/prediction and will be call Internal
    my $stop = 0;
    my $start = 0;
    my $numstops = 0;
    my $numstarts = 0;
    foreach my $gffrec (@{$transcript}){
	$numexon++ if $gffrec->[2] eq "exon";
	$numcds++ if $gffrec->[2] eq "CDS";
	if ($gffrec->[2] eq "start_codon"){
	  #print STDERR "start $gffrec->[2]\n";
	  $numstarts++; 
	  if ($gffrec->[6] eq '+'){
	    $start = $gffrec->[3];
	  }else{
	    $start = $gffrec->[4];
	  }
	}
	if ($gffrec->[2] eq "stop_codon"){
	  $numstops++; 
	  if ($gffrec->[6] eq '+'){
	    $stop = $gffrec->[4];
	  }else{
	    $stop = $gffrec->[3];
	  }
	}
    }
    
    foreach my $gffrec (@{$transcript}){
      $gffrec->[8]=~s/;$//;
	if ($gffrec->[3] < $tgff[3]){
	  $tgff[3]=$gffrec->[3];
				     #print STDERR "Changed tx coord\n";
	}
	if ($gffrec->[4] > $tgff[4]){
	  $tgff[4]=$gffrec->[4];#print STDERR "Changed tx coord\n";
	}
	if ($gffrec->[2] eq "CDS"){
	  if($ov == 2){
	    $gffrec->[8]="$tid";
	  }else{
	    $gffrec->[8].=";ID=$tid"."_cds";
	  }
	    push @cds, $gffrec;
	}elsif($gffrec->[2] eq "exon"){
	  if($ov == 2){
	    $gffrec->[8]="$tid";
	  }else{
	    if ($gffrec->[6] eq "+"){
		$gffrec->[8].=";ID=$tid"."_exon$count";
	    }else{
		$gffrec->[8].=";ID=$tid"."_exon".($numexon-$count + 1);
	    }
	    $count++;
	  }
	}else{
	  $gffrec->[7]='.';
	}
    }
    if($ov != 2){
      if (exists $tx{$tid}){
	print (join("\t",@{$tx{$tid}}),"\n");
      }else{
	print (join("\t",@tgff),"\n");
      }
    }
    my @scds = sort gffsort @cds; #make sure cds records are sorted
      
    
    ### Check for stop

    if ($check_start_stop_codon) {
      next if $fl && !($numstarts && $numstops);
      if ($addstop && @scds && $numstops == 1) {
	if ($tgff[6] eq "+") {
	  $scds[-1]->[4]+=3;
	} else {
	  $scds[0]->[3]-=3;
	}
      }
      if ($fits) {
	if ($numcds == 1) {
	  if ($numstarts && $numstops) {
	    $scds[0]->[2] = 'Single';
	  } elsif ($numstarts) {
	    $scds[0]->[2] = 'First';
	  } elsif ($numstops) {
	    $scds[0]->[2] = 'Terminal';
	  } else {
	    $scds[0]->[2] = 'Internal';
	  }
	} else {
	  for (my $i = 0;$i<$numcds;$i++) {
	    $scds[$i]->[2] = 'Internal'; 
	    if (($i == 0)&&($scds[$i]->[6] eq '+')) {
	      $scds[$i]->[2] = 'First' if ($numstarts && $scds[$i]->[3] == $start);
	    } else {
	      $scds[$i]->[2] = 'Terminal' if ($numstops && $scds[$i]->[3] == $stop);
	    }
	    if (($i ==  ($numcds -1))&&($scds[$i]->[6] eq '+')) {
	      $scds[$i]->[2] = 'Terminal' if ($numstops && $scds[$i]->[4] == $stop);
	    } else {
	      $scds[$i]->[2] = 'First' if ($numstarts && $scds[$i]->[4] == $start);
	    }
	  }
	} 
      }
    } else {
      if ($addstop && @scds) {
	if ($tgff[6] eq "+") {
	  $scds[-1]->[4]+=3;
	} else {
	  $scds[0]->[3]-=3;
	}
      }
      if ($fits) {
	if ($numcds == 1) {
	  $scds[0]->[2] = 'Single';
	} else {
	  for (my $i = 0;$i<$numcds;$i++) {
	    $scds[$i]->[2] = 'Internal';
	    $scds[$i]->[2] = 'First' if $i == 0;
	    $scds[$i]->[2] = 'Terminal' if $i == ($numcds -1);
	  }

	}
      
      }
    }
    if ($removestop && @scds) {
      if ($tgff[6] eq "+") {
	$scds[-1]->[4]-=3;
      } else {
	$scds[0]->[3]+=3;
      }
    }
    foreach my $gffrec (@{$transcript}) {
      if ($ov != 2) {
	$gffrec->[8]=~s/;$//;
	print (join("\t",@$gffrec),"\n");
       } else {
       	if ($gffrec->[2] =~ m/exon|First|Internal|Terminal|Single|CDS/) {
       	  print (join("\t",@$gffrec),"\n");
       	}
       }
    }
    print "###\n" if ($ov!=2);
  }

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
