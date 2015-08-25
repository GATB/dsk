# DSK  short manual

## Compilation

If you retrieved a source archive, you can use cmake to compile the tool:

    mkdir build; cd build; cmake ..; make 


## Usage

* type `./dsk` for usage.


## Input

File input format can be fasta, fastq, either gzipped or not.

To pass several files as input, separate file names by a comma (","), for example:  

    ./dsk  -file A1.fa,A2.fa,A3.fa  -kmer-size 31

Alternatively, input can be a list of files (one file per line): 

    ls -1 *.fastq > list_reads
    ./dsk -file list_reads
    

## Results visualisation

DSK uses a binary format for its output. 

By default, the format is HDF5, so you can use HDF5 tools for getting ascii output
from the HDF5 output (such tools are provided with DSK distribution).

In particular, the HDF5 output of DSK holds two sets of data:
   
* the couples of (kmer,abundance)
* the histogram of abundances
    
You can get content of a dataset (kmers in partition 0, and whole histogram) with:
    
    h5dump -y -d dsk/solid/0          output.h5
    h5dump -y -d histogram/histogram  output.h5

To see the results as a list of "[kmer] [count]\n", use the `dsk2ascii` binary with the output file from dsk 

    dsk2ascii -file output.h5 -out output.txt
    
To plot kmer coverage distribution,    
    
    h5dump -y -d histogram/histogram  output.h5  | grep "^\ *[0-9]" | tr -d " " | paste - - | gnuplot -p -e 'plot  "-" with lines'     


## Kmers and their reverse compelments

DSK converts kmers to a canonical representation with respect to reverse-complementation. 

This means that a kmer and its reverse complement are considered to be the same object. For example, with k=3 and assuming the kmer AAA and its reverse complement TTT are both present the input dataset, then DSK will consider that one of them is the canonical kmer, e.g. AAA. If AAA is present 2 times and TTT 3 times, then DSK will return that the count of AAA is 5 and wont return the count of TTT at all.

Note: the canonical kmer is not the lexicographically smaller one! DSK uses a different ordering for faster performance. Specifically, DSK considers tha A<C<T<G and returns the lexicographically smaller kmer using this alphabet order. So, in the example above, AAA is indeed the canonical kmer. For the GTA/TAC pair, the lexicographically smallest is GTA however the canonical kmer is TAC (as DSK considers that T<G).


## Larger k-mer sizes

DSK supports arbitrary large k-mer lengths. 
Just compile from the source, to support k-mer lengths up to, say, 160, type this in the build folder:

    rm -Rf CMake* && cmake -Dk4=160 .. && make

