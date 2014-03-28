# DSK  short manual

## Usage

* type `./dsk` for usage.


## Input

File input can be fasta, fastq, gzipped or not.

To pass several files as input : pass a list of file names (separator is ,) ex:  

    ./dsk  -file A1.fa,A2.fa,A3.fa  -kmer-size 31
    

## Results visualisation

DSK uses a binary format for its output. 

By default, the format is HDF5, so you can use HDF5 tools for getting ascii output
from the HDF5 output (such tools are provided with DSK distribution).

In particular, the HDF5 output of DSK holds two sets of data:
   
* the couples of (kmer,abundance)
* the histogram of abundances
    
You can get content of a dataset (description + data) with:
    
    h5dump -y -d dsk/solid      output.h5
    h5dump -y -d dsk/histogram  output.h5

To see the results as a list of "[kmer], [count]\n", type 

    h5dump -y -d dsk/solid output.h5  | grep "^\ *[0-9]" | tr -d " " | paste - -
    
To plot kmer coverage distribution,    
    
    h5dump -y -d dsk/histogram  output.h5  | grep "^\ *[0-9]" | tr -d " " | paste - - | gnuplot -p -e 'plot  "-" with lines'     

    