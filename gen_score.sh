#!/bin/bash

if [ "$#" -ne 3 ]
then
echo "Usage: ./gen_score.sh <probname> <score_type> <bf>
where <probname> is the name of the csv file in ./CSV/,
<score_type> is either \"BIC\" or \"BDeu\",
and <bf> is the desired Bayes factor,
e.g., ./gen_score.sh wine BIC 3 for ./CSV/wine.csv, BIC scores and a Bayes factor of 3.
"
exit 1
fi

mkdir -p ./scores/settings

echo "#GOBNILP parameters for generating scores
gobnilp/delimiter = \",\"
gobnilp/mergedelimiters = FALSE
gobnilp/outputfile/scores = \"./scores/$1.$2.$3\"
gobnilp/scoring/arities = FALSE
gobnilp/scoring/names = FALSE
gobnilp/scoring/prune = TRUE
gobnilp/scoring/prunegap = -$3
gobnilp/scoring/score_type = \"$2\"" > ./scores/settings/$1.$2.$3

if [ $2 == "BIC" ]
then
num_lines=$(wc -l < ./CSV/$1.csv)
num_data=$(($num_lines-2))
num_par=$(echo "l($num_data)/l(2)+l($3)" | bc -l)
num_par=$(echo "($num_par+1)/1" | bc)

echo "gobnilp/scoring/palim = $num_par" >> ./scores/settings/$1.$2.$3
fi

./gobnilp/bin/gobnilp -x -g=./scores/settings/$1.$2.$3 -f=dat ./CSV/$1.csv
