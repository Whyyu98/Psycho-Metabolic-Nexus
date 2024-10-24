#!/bin/bash

generate_and_run_cmd() {
  line=$1
  line=$(echo "$line" | xargs)  
  
  group1=$(echo "$line" | cut -d'-' -f1)
  group2=$(echo "$line" | cut -d'-' -f2)
  
  cmd="./mtag.py --sumstats ./PD/${group1}.txt,./MeTs/${group2}.txt --out ./PD-MeTs/${group1}-${group2} --stream_stdout --perfect_gencov --equal_h2 --force"
  
  echo ": $cmd"
  eval $cmd
}

export -f generate_and_run_cmd

cat Group.txt | tr -d '\r' | parallel -j 20 generate_and_run_cmd
