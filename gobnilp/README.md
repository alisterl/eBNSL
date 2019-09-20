(For information on the Python implementation of GOBNILP please consult the README.md file
in the Python directory. This readme only concerns the C implementation.)

GOBNILP is a C program which uses SCIP ( http://scip.zib.de/ ) to
learn Bayesian network structure from complete discrete data (or
precomputed local scores). If the BLAS, LAPACK and LAPACKE libraries
are available it can also learn Gaussian BNs from continuous data (see
below for more details on that). The version of GOBNILP that exists in
this git repo is only guaranteed to work with SCIP 6.0.1. (It has yet
to be 'tuned' to account for the changes between SCIP 4.0.0 and SCIP
6.0.1, so performance is typically worse than was the case with SCIP
4.0.0.)  

When installing SCIP 6.0.1 install it using "make" rather than
"cmake". Let us suppose you have installed SCIP in the directory
~/local/src/scipoptsuite-6.0.1/scip. GOBNILP can either be installed
to learn from discrete data (or pre-computed local scores) only or to
additionally learn Gaussian networks using BGe scoring. In the latter
case it is necessary to have the BLAS, LAPACK and LAPACKE libraries
installed. Some advice about how to install these libraries is given
below, but first here is how to install GOBNILP to learn only from
discrete data. Just do:


```
#!bash

git clone https://bitbucket.org/jamescussens/gobnilp.git
cd gobnilp/
./configure.sh ~/local/src/scipoptsuite-6.0.1/scip
make
```

The 'configure.sh' script mentioned just above copies (and alters)
files from SCIP which cannot be included in this repo for licensing
reasons. If you have CPLEX installed (and accessible to SCIP) be sure
to do 'make LPS=cpx' instead of plain 'make' since this will lead to
considerably faster solving. 

To learn Gaussian networks GOBNILP uses BGe scoring to compute local
scores using the BLAS, LAPACK and LAPACKE libraries. We have found
installing OpenBLAS from source the easiest option (we have only
tried with Linux). Note that OpenBLAS includes LAPACK and LAPACKE
which is particularly convenient.

First go to https://www.openblas.net/ and download
the source. At time of writing this gets you the file
OpenBLAS-0.2.20.tar.gz. To install OpenBLAS successfully you will need a
Fortran compiler, such as gfortran, on your system. (If you are using Ubuntu
then "sudo apt-get install gfortran" should work.) To install the libraries in, say,
"/home/james/local/" do the following:

```
#!bash

tar zxvf OpenBLAS-0.2.20.tar.gz
cd OpenBLAS-0.2.20
make
make install PREFIX=/home/james/local
```

This will put the header file lapacke.h into /home/james/local/include
and libopenblas.a, the library itself, into
/home/james/local/lib. Obviously replace "/home/james/local/" with
whatever directory is convenient for you. Once you have done all this
you can install GOBNILP as follows:

```
#!bash

git clone https://bitbucket.org/jamescussens/gobnilp.git
cd gobnilp/
./configure.sh ~/local/src/scipoptsuite-6.0.0/scip /home/james/local
make
```

where again, of course, you can replace "/home/james/local/" with a
directory of your choice. Running configure.sh alters the GOBNILP
Makefile to work with OpenBLAS installed as above. If you are using a
Mac you may have problems running configure.sh since it uses 'sed' to
edit the Makefile which can fail on Macs. In this case just edit the
Makefile yourself to set the variable BLAS to TRUE and the variable
BLASDIR to your equivalent of "/home/james/local/", e.g. just add
two lines like this somewhere near the top of Makefile:

```
BLAS = TRUE
BLASDIR = /home/james/local
```


BGe scoring was implemented by Sam Vozza as a final-year project in
the Dept of Computer Science, University of York.


For full details on GOBNILP (including installation instructions for the stable release), please consult the GOBNILP page:

https://www.cs.york.ac.uk/aig/sw/gobnilp/

The manual for the stable release is here:

https://www.cs.york.ac.uk/aig/sw/gobnilp/manual.pdf

The manual for development version in this git repo can be found (as LaTeX source) in the manual directory.

If you wish to extend GOBNILP in some way it is best to exploit SCIP's modular architecture to do so. For example,

1. To allow the user to add a new sort of constraint to GOBNILP you should add a new "constraint handler" to GOBNILP. See e.g. cons_dagcluster.c in the source, this implements the key acyclicity constraint. Also see cons_ci.c which implements user-defined conditional independence constraints. For general details on constraint handlers in SCIP, see:
http://scip.zib.de/doc/html/CONS.php

2. To add a new heuristic algorithm (to generate better incumbent solutions ) you should add a "primal heuristic". GOBNILP currently has only one of these: heur_sinks.c. For general info on primal heuristics in SCIP see
http://scip.zib.de/doc/html/HEUR.php

ONGOING DEVELOPMENT WORK

1. Scott MacGregor has written an R package for GOBNILP which will be
distributed at some point.

TODO LIST

Here are some things which are missing from GOBNILP but which could be added:

1. Adding a pricer so that "family variables" are created only if they loosen ( ie increase ) the dual bound, rather than all at once at the beginning of solving.
