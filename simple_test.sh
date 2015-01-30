#!/bin/bash
#simple  test with synthetic data

################################################################################
echo -n "Testing single gz file ..........."
bin/dsk -file test/read50x_ref10K_e001.fasta.gz -kmer-size 27 -out test_dsk27 -verbose 0
bin/h5dump -y -d dsk/histogram test_dsk27.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_dsk27.histo

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
bin/dsk -file test/c1.fasta.gz,test/c2.fasta.gz,test/c3.fasta.gz,test/c4.fasta.gz  -kmer-size 27 -out test_dsk27 -verbose 0
bin/h5dump -y -d dsk/histogram test_dsk27.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_dsk27.histo

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
bin/dsk -file test/longread.fasta -kmer-size 27 -out test_long  -verbose 0
bin/h5dump -y -d dsk/histogram test_long.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_long.histo

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
bin/dsk -file test/shortread.fasta  -kmer-size 15  -abundance-min 1  -out test_short  -verbose 0
bin/dsk2ascii -file test_short -out test_short.parse_results  -verbose 0

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
bin/dsk -file test/shortread.fasta  -kmer-size 16 -out test_short  &> /dev/null

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
bin/dsk -file test/readN.fasta -kmer-size 20 -out test_N  -verbose 0
bin/h5dump -y -d dsk/histogram test_N.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > test_N.histo

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

