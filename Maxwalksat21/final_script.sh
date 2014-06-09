#!/bin/bash
infile=$1
outfile=$2
rm $outfile
for ((i=10;i<=200;i+=10))
  do
     for j in {1..5..1}
       do
         ./maxwalksat -maxtime $i -tries 5000000  < $infile | tail -5 | head -3 >> $outfile ;
         printf "\n" >> $outfile;
       done
 done
