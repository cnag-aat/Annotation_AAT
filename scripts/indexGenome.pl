#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use lib "/home/groups/assembly/talioto/myperlmods/";
use Bio::DB::Fasta;

my $fasta = 0;

GetOptions(
	   'f:s' => \$fasta
	  );
die "usage: $0 -f <fasta.fa>\n" if ! $fasta; 
my $db = Bio::DB::Fasta->new($fasta,'-reindex');

print STDERR "Finished indexing $fasta\n";
exit;
