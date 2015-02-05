# DSK  short manual

## Compilation

If you retrieved a source archive, you can use cmake to compile the tool:

    mkdir build; cd build; cmake ..; make 


## Usage

* type `./dsk` for usage.


## Input

File input can be fasta, fastq, gzipped or not.

To pass several files as input : pass a list of file names (separator is ,) ex:  

    ./dsk  -file A1.fa,A2.fa,A3.fa  -kmer-size 31

Alternatively, it takes a list of files (one file per line) as input: 

    ls -1 *.fastq > list_reads
    ./dsk -file list_reads
    

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

To see the results as a list of "[kmer] [count]\n", use the `dsk2ascii` binary with the output file from dsk 

    dsk2ascii -file output.h5 -out output.txt
    
To plot kmer coverage distribution,    
    
    h5dump -y -d dsk/histogram  output.h5  | grep "^\ *[0-9]" | tr -d " " | paste - - | gnuplot -p -e 'plot  "-" with lines'     


## Larger k-mer sizes

DSK supports arbitrary large k-mer lengths. 
Just compile from the source, to support k-mer lengths up to, say, 160, type this in the build folder:

    rm -Rf CMake* && cmake -Dk4=160 .. && make

