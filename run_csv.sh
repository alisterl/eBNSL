#!/bin/bash

if [ "$#" -ne 2 ]; then
echo "Usage: ./run_csv <probname> <bf>
where <probname> is the name of the csv file in ./CSV
and <bf> is the desired Bayes factor,
e.g., ./run_csv wine 3 for ./CSV/wine.csv and a Bayes factor of 3.
"
exit 1
fi

mkdir -p ./scores/settings
mkdir -p ./results/settings
sed -i "54s/.*/#define EXP_BF $2/" ./gobnilp/src/scoring.c

# Check the make command to add flags if it fails
make LPS=cpx -C ./gobnilp || exit 1

num_lines=$(wc -l < ./CSV/$1.csv)
num_data=$(($num_lines-2))
num_par=$(echo "l($num_data)/l(2)+l($2)" | bc -l)
num_par=$(echo "($num_par+1)/1" | bc)

echo "#GOBNILP parameters for generating scores
gobnilp/scoring/palim = $num_par
gobnilp/delimiter = \",\"
gobnilp/mergedelimiters = FALSE
gobnilp/scoring/names = FALSE
gobnilp/outputfile/scores = \"./scores/$1.$2\"
gobnilp/scoring/arities = FALSE" > ./scores/settings/$1.$2

./gobnilp/bin/gobnilp -x -g=./scores/settings/$1.$2 -f=dat ./CSV/$1.csv 

echo "#GOBNILP parameters for finding the optimal network
gobnilp/outputfile/solution = \"./results/$1.$2.opt\"" > ./results/settings/$1.$2.opt

./gobnilp/bin/gobnilp -f=jkl -g=./results/settings/$1.$2.opt ./scores/$1.$2

score=$(tail -1 ./results/$1.$2.opt | awk '{print $4}')
score=$(echo "$score-l($2)" | bc -l)

echo "#GOBNILP parameters for collecting networks
gobnilp/objlimit = $score
gobnilp/outputfile/countsols = \"./results/$1.$2\"
gobnilp/countsols = TRUE
gobnilp/countsols/collect = TRUE
gobnilp/countsols/sollimit = 150000" > ./results/settings/$1.$2

./gobnilp/bin/gobnilp -g=./results/settings/$1.$2 -f=jkl ./scores/$1.$2
