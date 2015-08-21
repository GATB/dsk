#!/$bindir/bash
#simple  test with synthetic data

#look for the $bindir dir, normally it's in ../../$bindir because simple_test.sh is installed in share/dsk
#but it may be in ../$bindir if we're calling it from ./build
if [ -f "../../bin/dsk" ]
then
 bindir=../../bin
elif [ -f "./dsk" ]
then
 bindir=.
else
 echo "could not find a compiled dsk in ./$bindir nor ../../$bindir"
 exit 1
fi
echo "Testing $bindir/dsk"

################################################################################
echo -n "Testing single gz file ..........."
$bindir/dsk -file test/read50x_ref10K_e001.fasta.gz -kmer-size 27 -out test_dsk27 -max-memory 200 -verbose 0
$bindir/h5dump -y -d histogram/histogram test_dsk27.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_dsk27.histo

diff test_dsk27.histo  test/k27.histo > /dev/null

var=$?

if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
#    exit 1
fi

################################################################################
echo -n "Testing multiple gz files ........"
$bindir/dsk -file test/c1.fasta.gz,test/c2.fasta.gz,test/c3.fasta.gz,test/c4.fasta.gz -kmer-size 27 -out test_dsk27 -max-memory 200 -verbose 0
$bindir/h5dump -y -d histogram/histogram test_dsk27.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_dsk27.histo

diff test_dsk27.histo  test/k27.histo > /dev/null

var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
#    exit 1
fi

rm -f test_dsk27.*

################################################################################
echo -n "Testing long reads ..........."
$bindir/dsk -file test/longread.fasta -kmer-size 27 -out test_long  -verbose 0 -max-memory 200
$bindir/h5dump -y -d histogram/histogram test_long.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_long.histo

diff test_long.histo  test/rlong.histo > /dev/null

var=$?

if [ $var -eq 0 ]
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
#    exit 1
fi
rm -f test_long.*


################################################################################
echo -n "Testing k = readlen ..........."
$bindir/dsk -file test/shortread.fasta  -kmer-size 15  -abundance-min 1  -out test_short  -verbose 0 -max-memory 200
$bindir/dsk2ascii -file test_short -out test_short.parse_results  -verbose 0

diff test_short.parse_results  test/short.parse_results > /dev/null

var=$?

if [ $var -eq 0 ]
then
    echo  PASSED
else
    echo  FAILED
fi

rm -f test_short.*

################################################################################
echo -n "Testing k = readlen+1 ..........."
$bindir/dsk -file test/shortread.fasta  -kmer-size 16 -out test_short -max-memory 200  &> /dev/null

[ -s test_short.parse_results ]
var=$?

if [ $var -eq 1 ]
then
    echo  PASSED
else
    echo  FAILED
fi

rm -f test_short.*

################################################################################
echo -n "Testing read with N ........"
$bindir/dsk -file test/readN.fasta -kmer-size 20 -out test_N  -verbose 0 -max-memory 200
$bindir/h5dump -y -d histogram/histogram test_N.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_N.histo

diff test_N.histo  test/readN.histo > /dev/null

var=$?

if [ $var -eq 0 ]
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
#    exit 1
fi

rm -f test_N.*

