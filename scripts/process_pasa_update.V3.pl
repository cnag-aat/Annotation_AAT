#!/usr/bin/perl

############################################################
#
# script to process the pasa uodates comparing them with the EVM models in order to add some models or changing the updates in some cases.
#
# Author: J. Gómez Garrido
############################################################

use strict;
use warnings;

my $update = $ARGV[0];
defined($update) || die ("##ERROR## This script requires a pasa update file as first argument\n");
my $evm_models = $ARGV[1];
defined($evm_models) || die ("##ERROR## This script requires an evm output file as third argument.\n");

open UPDATE, "<", "$update";

my %genes; 
my %transcripts;
my %exons;
my %cds;
my $cds;
my $gene;
my $mRNA;
my $mRNA_parent;
my $cds_parent;
my $len;
my %updt_coords;
my $n;
while (<UPDATE>) {
    chomp;
    next if /^\#/o;
    next if /^$/o;
    #print "$_\n";
    $_.=";";
    my @line = split /\t/, $_;
    
    if ($line[2] eq 'gene'){
        if ($line[8] =~ m/ID=([^;]+)/) {$gene = $1;}
      #  print "$gene\n";
        $genes{$gene}=$_;
    #    print "$gene\n";
    }
    elsif ($line[2] eq 'mRNA'){
        if ($line[8] =~ m/ID=([^;]+)/) {$mRNA = $1;}
        if ($line[8] =~ m/Parent=([^;]+)/) {$mRNA_parent = $1;}
     #   print "$mRNA\t$mRNA_parent\n";
        $transcripts{$mRNA_parent}{$mRNA}=$_;
        $n = 0;
      #  print MRNAS "$_\n";
    }
    else {
        if ($line[8] =~ m/Parent=([^;]+)/) {$cds_parent = $1;}
        if ($line[8] =~ m/ID=([^;]+)/) {$cds = $1;}
      #  print "$cds_parent\n";
        $exons{$cds_parent} = "" if (!exists $exons{$cds_parent}); 
        $exons{$cds_parent} = $exons{$cds_parent} . $_ . "\n";
        if ($line[2] eq 'CDS') {
            $len=$line[4]-$line[3];
            $cds{$gene}{$cds_parent}->{lengths} = 0 if (!exists $cds{$gene}{$cds_parent}->{lengths});
            $cds{$gene}{$cds_parent}->{lengths}= $cds{$gene}{$cds_parent}->{lengths} + $len;
            $updt_coords{$gene}{$cds_parent}{$n}->{start}=$line[3];
            $updt_coords{$gene}{$cds_parent}{$n}->{end}=$line[4];
            $n++;
           # print  "$gene\t$cds_parent\t$cds\t$len\t$cds{$gene}{$cds_parent}->{lengths}\n";
        }
    }
}
close UPDATE;

open EVM, "<", "$evm_models"; 
my $evm_gene;
my $cds_model;
my %mrnas;
my %evm_genes;
my $evm_mRNA;
my %evm_transcripts;
my $evm_cds;
my $evm_cds_parent;
my %evm_exons;
my %evm_cds;
my %evm_coords;
my %evm_mrna_coords;
while (<EVM>){
    chomp;
    next if /^\#/o;
    next if /^$/o;
    #print "$_\n";
    $_ .= ";";
    my @line = split /\t/, $_;
    if ($line[2] eq 'gene'){
        if ($line[8] =~ m/ID=([^;]+)/) {
            $evm_gene = $1;
            $evm_genes{$evm_gene}=$_;
           # print STDERR "$_\n";
        }
    #print "$evm_gene\n";
    }
    elsif ($line[2] eq 'mRNA'){
        if ($line[8] =~ m/ID=([^;]+)/) {$evm_mRNA = $1;}
        if ($line[8] =~ m/Parent=([^;]+)/) {
            $evm_gene = $1;            
        }
        elsif ($line[8] =~ m/Parent=([^;]+)/) {
            $evm_gene = $1;            
        }
        $evm_transcripts{$evm_gene}{$evm_mRNA} = $_; 
        $n = 0;
        $evm_mrna_coords{$evm_gene}{$evm_mRNA}->{start} = $line[3];
        $evm_mrna_coords{$evm_gene}{$evm_mRNA}->{end} = $line[4];
        #print STDERR "$evm_gene\t$evm_mRNA\n";
    }
    else {
        if ($line[8] =~ m/Parent=([^;]+)/) {$evm_cds_parent = $1;}
        if ($line[8] =~ m/ID=([^;]+)/) {$evm_cds = $1;}
      #  print "$cds_parent\n";
        $evm_exons{$evm_cds_parent} = "" if (!exists $evm_exons{$evm_cds_parent}); 
        $evm_exons{$evm_cds_parent} = $evm_exons{$evm_cds_parent} . $_ . "\n";
        if ($line[2] eq 'CDS') {
            $len=$line[4]-$line[3];
            $evm_cds{$evm_gene}{$evm_cds_parent}->{lengths} = 0 if (!exists $evm_cds{$evm_gene}{$evm_cds_parent}->{lengths});
            $evm_cds{$evm_gene}{$evm_cds_parent}->{lengths}= $evm_cds{$evm_gene}{$evm_cds_parent}->{lengths} + $len;
            $evm_coords{$evm_gene}{$evm_cds_parent}{$n}->{start}=$line[3];
            $evm_coords{$evm_gene}{$evm_cds_parent}{$n}->{end}=$line[4];
            $n++;
           # print  "$evm_gene\t$evm_cds_parent\t$evm_cds\t$len\t$evm_cds{$evm_gene}{$evm_cds_parent}->{lengths}\n";
        }
    }
}
close EVM;
my %final_mrna;
my %evm_length;
my %updt_length;
my %final_mrna_comb;
my %print_evm_comb;
my %final_gene_comb;
my %conserved_cds;
my %print_mrna_evm;
my %ids;
my $tr;
foreach my $parse_gene (keys %genes) {
  #  my $size = scalar keys %{$evm_cds{$parse_gene}};
  #  print STDERR "$parse_gene\t$size\n";
  # print STDERR "$parse_gene\n";
   if (!exists $evm_cds{$parse_gene}) {
       my @test = split /_evm/, $parse_gene;
       my $total = scalar @test;
       my $i = 0;
    #   print STDERR "$parse_gene\n";
       $evm_length{$parse_gene} = 0;
    #  print STDERR "$total\n";
       while ($i < $total) {
           if ($i != 0) {
               $tr = "evm" . $test[$i]; 
           }
           else {
               $tr = $test[$i];
           }
           $i = $i + 1;
     #      print STDERR "$tr\n";
           foreach (keys %{$evm_transcripts{$tr}}) {
               #print STDERR "$evm_transcripts{$tr}{$_}\n";
               $evm_length{$parse_gene} = $evm_length{$parse_gene} + $evm_cds{$tr}{$_}->{lengths};
              # print STDERR "$tr\t$_\t$parse_gene\t$evm_length{$parse_gene}\n";
           }
      }
      foreach my $updt_mrna (keys %{$cds{$parse_gene}}) {
         #print STDERR "$parse_gene\t$evm_length{$parse_gene}\t$cds{$parse_gene}{$updt_mrna}->{lengths}\n";
         my $cov = $cds{$parse_gene}{$updt_mrna}->{lengths} / $evm_length{$parse_gene};
         if ($cov > 0.7) {
             $final_mrna_comb{$parse_gene}{$updt_mrna}++;
             my $i = 0;
             my @tr;
             while ($i < $total) {
              # my @tr = ("$test[$i]", "$test[$i+1]");
                push (@tr, "$test[$i]" );
               # print STDERR "@tr\n";
                $i++;
            }
            my $tr = join "_evm", @tr;
            #    print STDERR "$tr\n";
            #   $i = $i + 2;
             #  print STDERR "$evm_genes{$tr}\n";      
           foreach my $evm_con_cds (keys %{$evm_coords{$tr}}) {
               #  print STDERR "$evm_con_cds\n";
                   foreach my $evm_cds_checked(keys %{$evm_coords{$tr}{$evm_con_cds}}) {
                   #  print STDERR "$parse_gene\t$evm_con_cds\t$evm_cds_checked\n";
                      my $start = $evm_coords{$tr}{$evm_con_cds}{$evm_cds_checked}->{start};
                      my $end = $evm_coords{$tr}{$evm_con_cds}{$evm_cds_checked}->{end};
                      foreach my $updt_con_cds(keys %{$updt_coords{$parse_gene}}) {
                      #     print STDERR "$updt_con_cds\n";
                          foreach (keys %{$updt_coords{$parse_gene}{$updt_con_cds}}) {
                            # print STDERR "{$evm_con_cds}\t$evm_cds_checked\t$updt_con_cds\t$_\n";
                             if ($updt_coords{$parse_gene}{$updt_con_cds}{$_}->{start} == $start && $updt_coords{$parse_gene}{$updt_con_cds}{$_}->{end} == $end) {
                                  $conserved_cds{$evm_con_cds}{$evm_cds_checked}++;
                                  #print STDERR "{$evm_con_cds}\t$evm_cds_checked\n";
                             }
                          }              
                          #$print_mrna_evm{$evm_con_cds}++ if (!exists $print_mrna_evm{$evm_con_cds} && !exists $conserved_cds{$evm_con_cds}{$evm_cds_checked});
                        #  $ids{$evm_con_cds}++ if (!exists $conserved_cds{$evm_con_cds}{$evm_cds_checked});
                        #  print STDERR "$evm_con_cds}\n" if (!exists $conserved_cds{$evm_con_cds}{$evm_cds_checked});
                      }
                   #  print STDERR "$evm_con_cds}\n" if (!exists $conserved_cds{$evm_con_cds}{$evm_cds_checked});
                      $ids{$tr}++ if (!exists $ids{$tr} && !exists $conserved_cds{$evm_con_cds}{$evm_cds_checked});
                   }
               }              
           #  }
             # print STDERR "$updt_mrna\n";
         }
         else {
            # print STDERR "@test\n";
            my $i = 0;
            while ($i < $total) {
            #   my @tr = ("$test[$i]", "$test[$i+1]");
               # print STDERR "@tr\n";
             #  my $tr = join "_", @tr;
               # print STDERR "$tr\n";
               my $tr = $test[$i];
               $i = $i + 1;
               #print STDERR "$evm_genes{$tr}\n";
               $print_evm_comb{$tr}++ if (!exists $print_evm_comb{$tr});
            }        
         }  
     }
   }    
   foreach my $updt_mrna (keys %{$cds{$parse_gene}}) {
      #print STDERR "$updt_mrna\n";
      #print STDERR "$updt_mrna\n" if (!exists $evm_cds{$parse_gene});
      foreach my $evm_mrna (keys %{$evm_cds{$parse_gene}}) {
         my $updt_length = $cds{$parse_gene}{$updt_mrna}->{lengths};
         my $evm_length = $evm_cds{$parse_gene}{$evm_mrna}->{lengths};
       #  print "$parse_gene\t$updt_mrna\t$evm_mrna\t$updt_length\t$evm_length\n";
         my $cov = $updt_length / $evm_length;
         if ($cov > 0.7) {
            $final_mrna{$parse_gene}{$updt_mrna}++;
            foreach my $evm_con_cds (keys %{$evm_coords{$parse_gene}}) {
             #  print STDERR "$evm_con_cds\n";
               foreach my $evm_cds_checked(keys %{$evm_coords{$parse_gene}{$evm_con_cds}}) {
                #  print STDERR "$parse_gene\t$evm_con_cds\t$evm_cds_checked\n";
                  my $start = $evm_coords{$parse_gene}{$evm_con_cds}{$evm_cds_checked}->{start};
                  my $end = $evm_coords{$parse_gene}{$evm_con_cds}{$evm_cds_checked}->{end};
                  foreach my $updt_con_cds(keys %{$updt_coords{$parse_gene}}) {
                  #     print STDERR "$updt_con_cds\n";
                      foreach (keys %{$updt_coords{$parse_gene}{$updt_con_cds}}) {
                         #  print STDERR "{$evm_con_cds}\t$evm_cds_checked\t$updt_con_cds\t$_\n";
                           if ($updt_coords{$parse_gene}{$updt_con_cds}{$_}->{start} == $start && $updt_coords{$parse_gene}{$updt_con_cds}{$_}->{end} == $end) {
                               $conserved_cds{$evm_con_cds}{$evm_cds_checked}++;
                             #  print STDERR "{$evm_con_cds}\t$evm_cds_checked\n";
                           }              
                       }

                    }
                    $print_mrna_evm{$evm_con_cds}++ if (!exists $print_mrna_evm{$evm_con_cds} && !exists $conserved_cds{$evm_con_cds}{$evm_cds_checked});

                }              
            }
          #  print "$parse_gene\t$updt_mrna\t$evm_mrna\t$cov\n";
         }
      }
   }
}

foreach (keys %print_mrna_evm) {
 # print STDERR "$_\n";
}

foreach my $print (keys %genes){
    if (exists $final_mrna{$print}) {
        my @line = split /\s+/, $genes{$print};
        my $start = $line[3];
        my $end = $line[4];
        my $change_start = 0;
        my $change_end = 0;
        my $mrna_evm_print;
        foreach (keys %{$evm_transcripts{$print}}) {
            $mrna_evm_print = $_;
            if (exists $print_mrna_evm{$mrna_evm_print}) {
                if ($evm_mrna_coords{$print}{$mrna_evm_print}->{start} < $start) {
                   $change_start = 1;
                }
                if ($evm_mrna_coords{$print}{$mrna_evm_print}->{end} > $end) {
                   $change_end = 1;
                }
            }
        }
        if ($change_start == 1 || $change_end == 1) {
            my $print_start;
            my $print_end;
            if ($change_start == 1) {
                $print_start = $evm_mrna_coords{$print}{$mrna_evm_print}->{start};
            }
            else {
                $print_start = $start;
           }
           if ($change_end == 1) {
               $print_end = $evm_mrna_coords{$print}{$mrna_evm_print}->{end};
           }
           else {
               $print_end = $end;
           }
           print "$line[0]\t$line[1]\t$line[2]\t$print_start\t$print_end\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
        }
        else {
           print "$genes{$print}\n"; 
        }
        foreach (keys %{$transcripts{$print}}) {
            print "$transcripts{$print}{$_}\n" if (exists $final_mrna{$print}{$_});
            print "$exons{$_}" if (exists $final_mrna{$print}{$_});
        }
        foreach (keys %{$evm_transcripts{$print}}) {
            my @line = split /\t/, $evm_transcripts{$print}{$_};
            if ($line[8] =~ m/ID=([^;]+)/) {$mRNA = $1;}
          #  print "$mRNA\n";
            $evm_transcripts{$print}{$_} =~ s/;Parent=/.evm;Parent=/;
          #  print "$evm_exons{$_}\n" if (exists $print_mrna_evm{$_});
            $evm_exons{$_} =~ s/Parent=([^;]+)/Parent=$1.evm/g;
          #  print "$evm_exons{$_}\n" if (exists $print_mrna_evm{$_});
            $evm_exons{$_} =~ s/ID=cds.([^;]+)/ID=cds.$1.evm/g;
            print "$evm_transcripts{$print}{$_}\n" if (exists $print_mrna_evm{$_});
            print "$evm_exons{$_}" if (exists $print_mrna_evm{$_});
        }
    }
    else {
       $ids{$print}++;
    }
}

foreach my $test (keys %final_mrna_comb) {
   print "$genes{$test}\n";
   foreach (keys %{$final_mrna_comb{$test}}) {
      print "$transcripts{$test}{$_}\n";
      print "$exons{$_}";
   }
  # print STDERR "$test\n";
}

foreach my $print_evm (keys %evm_genes) {
    if (exists $ids{$print_evm} || (exists $print_evm_comb{$print_evm} && !exists $final_gene_comb{$print_evm})) {
  #    #  print STDERR "$print_evm\n";
        print "$evm_genes{$print_evm}\n";
        foreach (keys %{$evm_transcripts{$print_evm}}) {
           print "$evm_transcripts{$print_evm}{$_}\n";
           print "$evm_exons{$_}";
        }
    } 
}

