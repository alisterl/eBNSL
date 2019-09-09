#!/bin/bash

if [ "$#" -ne 2 ]; then
echo "Usage: ./run_score <scorename> <bf>
where <scorename> is the name of the scoring file in ./scores
and <bf> is the desired Bayes factor,
e.g., ./run_csv wine.3 3 for ./scores/wine.3 and a Bayes factor of 3.
"
exit 1
fi

mkdir -p ./results/settings

echo "#GOBNILP parameters for finding the optimal network
gobnilp/outputfile/solution = \"./results/$1.$2.opt\"" > ./results/settings/$1.$2.opt

./gobnilp/bin/gobnilp -f=jkl -g=./results/settings/$1.$2.opt ./scores/$1

score=$(tail -1 ./results/$1.$2.opt | awk '{print $4}')
score=$(echo "$score-l($2)" | bc -l)

echo "#GOBNILP parameters for collecting networks
gobnilp/objlimit = $score
gobnilp/outputfile/countsols = \"./results/$1.$2\"
gobnilp/countsols = TRUE
gobnilp/countsols/collect = TRUE
gobnilp/countsols/sollimit = 150000" > ./results/settings/$1.$2

./gobnilp/bin/gobnilp -g=./results/settings/$1.$2 -f=jkl ./scores/$1
