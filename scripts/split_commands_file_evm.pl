#!/usr/bin/perl

use strict;
use warnings;

############################################################################################################################################################
## This script splits the command file written by evm in 16 files in order to the commands be run in 16 CPUs. The script evm.sh uses it, but if you want to ## call it you must write:
##       ./split_commands_file_evm.pl evm.cmd
##
## Author: Jèssica Gómez
## Date: 13032014
##
###########################################################################################################################################################

my $cmdfile = $ARGV[0];
defined($cmdfile) || die ("##ERROR## This script requires a command file as argument\n");
my $fragments = $ARGV[1];
open INFILE, "<", "$cmdfile" || die ("##ERROR## Cannot open file $cmdfile\n");

my @names = split /\./, $cmdfile;
my $name = $names[0];

my $i = 0;
my $j = 0;
my $c = 1;
my @file = <INFILE>;
my $lv = scalar(@file);


open OUTFILE, ">", "$name.$c.cmd";


while ($i < $fragments) { 
   while ($j < $lv) {
      my @line = split /\s+/, $file[$j];
      chomp @line;
      select OUTFILE;
      print "@line\n";
      $j = $j + $fragments;
    }
    $i++;
    $j = $i;
    close OUTFILE;
    $c++;
    open OUTFILE, ">", "$name.$c.cmd" if $c <= $fragments;
}


