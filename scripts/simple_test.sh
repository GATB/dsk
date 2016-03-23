#!/bin/bash
#simple test with synthetic data

# look for dsk binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/dsk" ]
then
 bindir="../bin"
 h5bindir="../bin"
elif [ -f "../build/bin/dsk" ]
then
 bindir="../build/bin"
 h5bindir="../build/ext/gatb-core/bin"
else
 echo "could not find a compiled dsk binary"
 exit 1
fi

# look for test data directory. In devel mode, it's in ../test directory.
# In production mode, it's in . directory (relative to this script).
if [ -f "./read50x_ref10K_e001.fasta.gz" ]
then
 testdir="."
elif [ -f "../test/read50x_ref10K_e001.fasta.gz" ]
then
 testdir="../test"
else
 echo "could not find test data directory"
 exit 1
fi

echo "Testing $bindir/dsk"

################################################################################
echo -n "Testing single gz file ..........."
$bindir/dsk -file $testdir/read50x_ref10K_e001.fasta.gz -kmer-size 27 -out test_dsk27 -max-memory 200 -verbose 0
$h5bindir/h5dump -y -d histogram/histogram test_dsk27.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_dsk27.histo

diff test_dsk27.histo  $testdir/k27.histo > /dev/null

var=$?
if [ $var -eq 0 ]
then
    echo  PASSED
else
    echo  FAILED
    exit 1
fi

################################################################################
echo -n "Testing multiple gz files ........"
$bindir/dsk -file $testdir/c1.fasta.gz,$testdir/c2.fasta.gz,$testdir/c3.fasta.gz,$testdir/c4.fasta.gz -kmer-size 27 -out test_dsk27 -max-memory 200 -verbose 0
$h5bindir/h5dump -y -d histogram/histogram test_dsk27.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_dsk27.histo

diff test_dsk27.histo  $testdir/k27.histo > /dev/null

var=$?
rm -f test_dsk27.*
if [ $var -eq 0 ]
then
    echo  PASSED
else
    echo  FAILED
    exit 1
fi


################################################################################
echo -n "Testing long reads ..........."
$bindir/dsk -file $testdir/longread.fasta -kmer-size 27 -out test_long  -verbose 0 -max-memory 200
$h5bindir/h5dump -y -d histogram/histogram test_long.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_long.histo

diff test_long.histo  $testdir/rlong.histo > /dev/null

var=$?
rm -f test_long.*
if [ $var -eq 0 ]
then
    echo  PASSED
else
    echo  FAILED
    exit 1
fi


################################################################################
echo -n "Testing k = readlen ..........."
$bindir/dsk -file $testdir/shortread.fasta  -kmer-size 15  -abundance-min 1  -out test_short  -verbose 0 -max-memory 200
$bindir/dsk2ascii -file test_short -out test_short.parse_results  -verbose 0

diff test_short.parse_results  $testdir/short.parse_results > /dev/null

var=$?
rm -f test_short.*
if [ $var -eq 0 ]
then
    echo  PASSED
else
    echo  FAILED
    exit 1
fi


################################################################################
echo -n "Testing k = readlen+1 ..........."
$bindir/dsk -file $testdir/shortread.fasta  -kmer-size 16 -out test_short -max-memory 200  &> /dev/null

[ -s test_short.parse_results ]
var=$?
rm -f test_short.*
if [ $var -eq 1 ]
then
    echo  PASSED
else
    echo  FAILED
    exit 1
fi


################################################################################
echo -n "Testing read with N ........"
$bindir/dsk -file $testdir/readN.fasta -kmer-size 20 -out test_N  -verbose 0 -max-memory 200
$h5bindir/h5dump -y -d histogram/histogram test_N.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_N.histo

diff test_N.histo  $testdir/readN.histo > /dev/null

var=$?
rm -f test_N.*
if [ $var -eq 0 ]
then
    echo  PASSED
else
    echo  FAILED
    exit 1
fi
