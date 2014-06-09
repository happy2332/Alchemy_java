infile='../student/student_aistats_29.cnf'
outfile='../student/student_aistats_29.results'
#rm $outfile
for i in {30..200..10}
  do
     for j in {1..5..1}
       do
         ./maxwalksat -maxtime $i -tries 5000000  < $infile | tail -5 | head -3 >> $outfile 
         printf "\n" >> $outfile
       done
 done
