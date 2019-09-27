#!/bin/bash

if [ "$#" -ne 2 ]
then
echo "Usage: ./run_score.sh <scorename> <bf>
where <scorename> is the name of the scoring file in ./scores/,
and <bf> is the desired Bayes factor,
e.g., ./run_score.sh wine.BIC.3 for ./scores/wine.BIC.3 and a Bayes factor of 3.
"
exit 1
fi

mkdir -p ./networks/settings

echo "#GOBNILP parameters for finding the optimal network
gobnilp/outputfile/solution = \"./networks/$1.opt\"" > ./networks/settings/$1.opt

./gobnilp/bin/gobnilp -f=jkl -g=./networks/settings/$1.opt ./scores/$1

score=$(tail -1 ./networks/$1.opt | awk '{print $4}')

echo "#GOBNILP parameters for collecting networks
gobnilp/countsols = TRUE
gobnilp/countsols/collect = TRUE
gobnilp/countsols/sollimit = 150000
gobnilp/objlimit = $(echo "$score-l($2)" | bc -l)
gobnilp/outputfile/countsols = \"./networks/$1\"" > ./networks/settings/$1

./gobnilp/bin/gobnilp -g=./networks/settings/$1 -f=jkl ./scores/$1
