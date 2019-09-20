# e - Bayesian Network Structure Learning

The package implements the following paper while the formal release of GOBNILP is under development. Please check the [GOBNILP website](https://www.cs.york.ac.uk/aig/sw/gobnilp/#aaai19) for updates.

Liao, Z. A., Sharma, C., Cussens, J., & van Beek, P. (2019, July). Finding all Bayesian network structures within a factor of optimal. In *Proceedings of the Thirty-Third AAAI Conference on Artificial Intelligence* (Vol. 33, pp. 7892-7899).
<details>
  <summary>BibTex</summary>  
<p>

```
@inproceedings{liao2019finding,
  title={Finding all {B}ayesian network structures within a factor of optimal},
  author={Liao, Zhenyu A. and Sharma, Charupriya and Cussens, James and van Beek, Peter},
  booktitle={Proceedings of the Thirty-Third AAAI Conference on Artificial Intelligence},
  volume={33},
  pages={7892--7899},
  year={2019}
}
```

</p>
</details>

## Getting Started

To use the package you need to first compile GOBNILP by following the steps provided in `./gobnilp/README.md`. You will also need to install [SCIP Optimization Suite](http://scip.zib.de/) required by GOBNILP.

### Prerequisites

* SCIP Optimization Suite 6.0.1 ([scipoptsuite-6.0.1](https://scip.zib.de/index.php#download))
* (Optional) IBM ILOG CPLEX Optimization Studio 12.9 ([free academic version](https://my15.digitalexperience.ibm.com/b73a5759-c6a6-4033-ab6b-d9d4f9a6d65b/dxsites/151914d1-03d2-48fe-97d9-d21166848e65/technology/data-science))

### Installing

1. (Optional) Install CPLEX
2. Install SCIP

  Type the following command in the `scipoptsuite-6.0.1` directory.
  * with CPLEX

    ```
    make LPS=cpx ZIMPL=false
    ```

    Please use the ***absolute path*** for CPLEX during the installation.
  * without CPLEX

    ```
    make ZIMPL=false
    ```

  If you encounter any problem, you can try adding some or all of the following flags.

  ```
  ZLIB=false GMP=false READLINE=false LPSOPT=opt-gccold OPT=opt-gccold
  ```

3. Compile GOBNILP

  The included GOBNILP is modified from the development version [[GitHash:db37374](https://bitbucket.org/jamescussens/gobnilp/src/db373747e76955f35437170f9641c9130bf50e9a/)]. Type the following command in `./gobnilp`.
  1. Link SCIP

    ```
    ./configure.sh <scipoptsuite-path>/scip
    ```

  2. Compile
    * with CPLEX

      ```
      make LPS=cpx ZIMPL=false
      ```

    * without CPLEX

      ```
      make ZIMPL=false
      ```
    You also need to add additional flags used in Step 2.

4. Check the installation

  You can use a small dataset to check whether GOBNILP has been properly installed. Type the following command in the package directory. It asks GOBNILP to find the optimal network using the pruned scoring file `wine.BIC.3`.

  ```
  ./gobnilp/bin/gobnilp ./scores/wine.BIC.3
  ```

  <details>
    <summary>Sample output</summary>  
  <p>

  ```
  GOBNILP version development [GitHash: 9f8daa2 ]
  Solving the BN structure learning problem using SCIP.

  SCIP version 6.0.2 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: CPLEX 12.9.0.0] [GitHash: e639a0059d]
  Copyright (C) 2002-2019 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)

  WARNING: Parameter file <gobnilp.set> not found - using default settings.
  WARNING: Input file format not recognised - assuming it is Jaakkola.
  File name:		./scores/wine.BIC.3
  Problem name:		wine
  Number of variables: 14
  Number of candidate parent sets: 592

  presolving (4 rounds: 4 fast, 3 medium, 2 exhaustive):
   99 deleted vars, 36 deleted constraints, 0 added constraints, 0 tightened bounds, 0 added holes, 81 changed sides, 81 changed coefficients
   0 implications, 4554 cliques
  presolved problem has 738 variables (738 bin, 0 int, 0 impl, 0 cont) and 450 constraints

  time | Best Network Found So Far |   Best Network Possible   | mem |  gap   |objleav|infleav
   0.2s|       -1.286816e+03       |        0.000000e+00       |9754k|    Inf |     0 |     0
   0.2s|       -1.276888e+03       |       -1.174130e+03       |9759k|   8.75%|     0 |     0
   0.2s|       -1.262882e+03       |       -1.194630e+03       |9767k|   5.71%|     0 |     0
   0.2s|       -1.262502e+03       |       -1.198862e+03       |9775k|   5.31%|     0 |     0
   0.2s|       -1.259041e+03       |       -1.214714e+03       |9820k|   3.65%|     0 |     0
   0.7s|       -1.258446e+03       |       -1.256731e+03       |  26M|   0.14%|     0 |     0
   1.4s|       -1.258446e+03       |       -1.258446e+03       |  50M|   0.00%|     0 |     0

  SCIP Status        : problem is solved [optimal solution found]
  Solving Time (sec) : 1.43
  Solving Nodes      : 1
  Primal Bound       : -1.25844590861000e+03 (13 solutions)
  Dual Bound         : -1.25844590861000e+03
  Gap                : 0.00 %
  0<-7,12, -121.535985
  1<-0, -72.056337
  2<-0, -90.136687
  3<-0, -123.583712
  4<-0,3, -94.653216
  5<-0, -116.096894
  6<-7, -54.213399
  7<-12, -71.276573
  8<-4,12, -93.494926
  9<-2,7, -95.559163
  10<-0, -62.365545
  11<-0, -94.357428
  12<- -123.759669
  13<-0, -45.356374
  BN score is -1258.445909

  ```

  </p>
  </details>

## Running the Experiments

* To generate Bayesian networks within a factor of optimal from datasets, type

  ```
  #!bash

  ./run_csv.sh <probname> <score_type> <bf>
  ```

  where `<probname>` is the name of the csv file in `./CSV/`, `<score_type>` is either `BIC` or `BDeu`, and `<bf>` is the desired Bayes factor. For example, type `./run_csv.sh wine BIC 3` for `./CSV/wine.csv`, `BIC` scoring function and a Bayes factor of `3`.

* To generate pruned scoring files from datasets, type

  ```
  #!bash

  ./gen_score.sh <probname> <score_type> <bf>
  ```

  where `<probname>` is the name of the csv file in `./CSV/`, `<score_type>` is either `BIC` or `BDeu`, and `<bf>` is the desired Bayes factor. For example, type `./gen_score.sh wine BDeu 3` for `./CSV/wine.csv`, `BDeu` scoring function and a Bayes factor of `3`.

* To collect Bayesian networks from pruned scoring files, type
  ```
  #!bash

  ./run_score <scorename> <bf>
  ```
  where `<scorename>` is the name of the scoring file in `./scores/` and `<bf>` is the desired Bayes factor. For example, type `./run_score.sh wine.BIC.3 3` for `./scores/wine.BIC.3` and a Bayes factor of `3`.

## Authors
Please send questions and comments to alister.liao AT uwaterloo.ca.
