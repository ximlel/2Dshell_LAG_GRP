#!/bin/bash  --login
#convert all .dat to .plt (tecplot data)
shopt -s expand_aliases
source ~/.bash_aliases

for i in $(ls *.dat)
do
preplot $i
done
