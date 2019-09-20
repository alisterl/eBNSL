#!/bin/bash

if [ "$#" -ne 3 ]
then
echo "Usage: ./run_csv.sh <probname> <score_type> <bf>
where <probname> is the name of the csv file in ./CSV/,
<score_type> is either \"BIC\" or \"BDeu\",
and <bf> is the desired Bayes factor,
e.g., ./run_csv.sh wine BIC 3 for ./CSV/wine.csv, BIC scores and a Bayes factor of 3.
"
exit 1
fi

./gen_score.sh $1 $2 $3

./run_score.sh $1.$2.$3 $3
