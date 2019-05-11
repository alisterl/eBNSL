# eBNSL

The package implements the following paper while the formal release is under development.
Please check GOBNILP website ( https://www.cs.york.ac.uk/aig/sw/gobnilp/#aaai19 ) for updates.

Liao, Z., Sharma, C., Cussens, J., & van Beek, P. (2018, December). Finding All Bayesian Network Structures within a Factor of Optimal. In The Thirty-Third AAAI Conference on Artificial Intelligence (AAAI-19). AAAI Press.

To use the package you need to first compile GOBNILP by following the steps provided in that folder. You will also need to install SCIP Optimization Suite ( http://scip.zib.de/ ) required by GOBNILP.

If you use any special flags when compiling GOBNILP, you need to add those flags to the make command in run_csv.sh.

To generate Bayesian networks from data, type
```
#!bash

./run_csv <probname> <bf>
```
where <probname> is the name of the csv file in ./CSV
and <bf> is the desired Bayes factor,
e.g., ./run_csv wine 3 for ./CSV/wine.csv and a Bayes factor of 3.

Please send questions and comments to alister.liao@uwaterloo.ca
