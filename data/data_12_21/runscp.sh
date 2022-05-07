#!/bin/bash 

for sigma in 1 5 10 15 20 25 30 33 35 37 40 45 46 47 48 49 50     #loop over parse variabe 2
do
  	scp `printf "mw936@amarel.rutgers.edu:/scratch/mw936/data_12_21/zips/sigma%03d.zip /Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/" $sigma`
done
