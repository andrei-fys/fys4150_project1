#!/bin/bash

#set -x

#  C++ programm should be run like this
#  ./gauss N computed calculated error thomas
#  ./LU N LU_file
#  where
#  N - number of grid points/matrix dimention
#  computed   - output file for Gauss elimination solution
#  calculated - output file for exact solution
#  error      - output file for relative error
#  thomas     - output file for Thomas algo. outtput
#  LU_file    - output file for LU algo. outtput

N=(10 100 1000 10000 100000 1000000)
#N=(10 100)
Gauss="Gauss_computed"
Thomas="Thomas_computed"
LU="LU_computed"
calculated="exact_solution"
error="relative_error"
results_dir="results"
time_file="time"

g++ -o gauss gauss.cpp
g++ -o LU LU.cpp lib.cpp

#make && make clean

for i in "${N[@]}";
do
	echo "Performing Gauss elimination for $i grid points"
	./gauss $i $Gauss"_"$i $calculated"_"$i $error"_"$i $Thomas"_"$i > $time_file"_"$i
	if [ ! -d $results_dir ]; then
		mkdir $results_dir
	fi
	if [ $i -le 1000 ]; then # LU for 10000 runs about 30 min on i7, so ...
		echo "Performing LU decomposition for $i grid points"
		./LU $i $LU"_"$i >> $time_file"_"$i
		mv $LU"_"$i $results_dir
	fi
	mv $Gauss"_"$i $calculated"_"$i $error"_"$i $Thomas"_"$i $time_file"_"$i $results_dir
done
