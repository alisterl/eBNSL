#! /bin/bash

function usage() {
   echo "USAGE"
   echo "Either: ./configure.sh <SCIPDIR>"
   echo "Or:     ./configure.sh <SCIPDIR> <BLASDIR>"
   echo "   where <SCIPDIR> is the path to the SCIP directory"
   echo "   where <BLASDIR>/include is the path to the BLAS include directory"
   echo "   and <BLASDIR>/lib is the path to the BLAS library directory"
   echo "   For example"
   echo "      ./configure.sh /opt/scipoptsuite-5.0.0/scip"
   echo "   or"
   echo "      ./configure.sh /home/userfs/j/jc33/local/src/scipoptsuite-5.0.1/scip/ /home/userfs/j/jc33/local/"
   exit -1
}

function errmsg() {
   echo $1
   exit -1
}

### Check there is the correct number of input arguments
if [ $# -ne 1 -a $# -ne 2 ]; then
   usage
fi

### Check the first input argument is a valid directory
SCIPDIR=$1
if [ ! -d $SCIPDIR ]; then
   errmsg "ERROR: Couldn't find \"$SCIPDIR\""
fi

if [ $# -eq 2 ]; then
    BLASDIR=$2
    BLASINCDIR="${BLASDIR}/include"
    BLASLIBDIR="${BLASDIR}/lib"
    if [ ! -d $BLASINCDIR ]; then
	errmsg "ERROR: Couldn't find BLAS include directory \"$BLASINCDIR\""
    fi
    if [ ! -d $BLASLIBDIR ]; then
	errmsg "ERROR: Couldn't find BLAS library directory \"$BLASLIBDIR\""
    fi
    SEDCMD="s|#OPEN BLAS not installed|#OPEN BLAS installed\nBLAS = TRUE\nBLASDIR = ${BLASDIR}\n|"
    echo "Editing Makefile"
    sed -i -e"${SEDCMD}" Makefile
fi

### Create the link to the SCIP directory
echo "Creating link to $SCIPDIR"
if [ -e scip ]; then
   errmsg "ERROR: Directory \"scip\" already exists"
fi
ln -s "$SCIPDIR" scip
if [ $? -ne 0 ]; then
   errmsg "ERROR: Couldn't create link to \"$SCIPDIR\""
fi

### Copy the linear ordering files across
echo "Copying linear ordering files"
if [ ! -e scip/examples/LOP/src/cons_lop.c ]; then
   errmsg "ERROR: Couldn't find file \"$SCIPDIR/examples/LOP/src/cons_lop.c\""
fi
cp scip/examples/LOP/src/cons_lop.c src/
if [ $? -ne 0 ]; then
   errmsg "ERROR: Couldn't copy \"$SCIPDIR/examples/LOP/src/cons_lop.c\""
fi
if [ ! -e scip/examples/LOP/src/cons_lop.h ]; then
   errmsg "ERROR: Couldn't find file \"$SCIPDIR/examples/LOP/src/cons_lop.h\""
fi
cp scip/examples/LOP/src/cons_lop.h src/
if [ $? -ne 0 ]; then
   errmsg "ERROR: couldn't copy \"$SCIPDIR/examples/LOP/src/cons_lop.h\""
fi

echo "Creating partial ordering files"
if [ ! -e partialc.patch ]; then
   errmsg "ERROR: Couldn't find patch file \"partialc.patch\""
fi
if [ ! -e partialh.patch ]; then
   errmsg "ERROR: Couldn't find patch file \"partialh.patch\""
fi

cp src/cons_lop.c src/cons_partialordering.c
if [ $? -ne 0 ]; then
   errmsg "ERROR: couldn't copy \"src/cons_lop.c\" to \"src/cons_partialordering.c\""
fi
patch src/cons_partialordering.c partialc.patch
if [ $? -ne 0 ]; then
   errmsg "ERROR: couldn't patch \"src/cons_partialordering.c\""
fi
cp src/cons_lop.h src/cons_partialordering.h 
if [ $? -ne 0 ]; then
   errmsg "ERROR: couldn't copy \"src/cons_lop.h\" to \"src/cons_partialordering.h\""
fi
patch src/cons_partialordering.h partialh.patch
if [ $? -ne 0 ]; then
   errmsg "ERROR: couldn't patch \"src/cons_partialordering.h\""
fi

### No problems occurred
echo "SUCCEEDED"
