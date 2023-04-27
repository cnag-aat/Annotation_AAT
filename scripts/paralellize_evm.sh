#!/bin/bash 

version=$1
split=$2
threads=$3
evm_script=$4

$split evm_$version.cmd $threads;
for ((i=1; i<=$threads; i++))
  do
  echo $i;
  $evm_script evm_$version.$i.cmd &
  done;

wait
