#!/usr/bin/env perl

my %gc; my $t=0;

while(<>){
  next if m/>/;
  chomp;
  
  while (/([ACGTacgt])/g){
    $gc{uc($1)}++;
    $t++;
  }
}
printf("%2.2f%%\n",($gc{G}/$t + $gc{C}/$t)*100);
