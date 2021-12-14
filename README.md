# DSK
| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/DSK-gitlab/job/tool-dsk-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/DSK-gitlab/job/tool-dsk-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/DSK-gitlab/job/tool-dsk-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/DSK-gitlab/job/tool-dsk-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is DSK? 
DSK is a k-mer counter for reads or genomes.

It takes as input a set of sequences in FASTA or FASTQ format (see "Input" section).
DSK outputs a set of solid kmers, i.e. kmers which occur more than a minimal amount of times in the input (-abundance-min parameter).
It also outputs the number of times these kmers occur.
See "Results visualization" section to learn how to use its output.

# Getting the latest source code

## Requirements

CMake 3.1+; see http://www.cmake.org/cmake/resources/software.html

C++/11 capable compiler (e.g. gcc 4.7+, clang 3.5+, Apple/clang 6.0+).

## Instructions

    # get a local copy of source code
    git clone --recursive https://github.com/GATB/dsk.git
    
    # compile the code an run a simple test on your computer
    cd dsk
    sh INSTALL

## Testing

If you retrieved the software from a git clone, do the following AFTER compiling it:

    cd ../scripts  # we suppose you are in the build directory
    ./simple_test.sh

If you retrieved the software from a binary distribution:

    cd test
    ./simple_test.sh

## Using

* type `./dsk` for usage instructions.


## Input

File input format can be fasta, fastq, either gzipped or not.

To pass several files as input, separate file names by a comma (","), for example:  

    ./dsk  -file A1.fa,A2.fa,A3.fa  -kmer-size 31

Alternatively, input can be a list of files (one file per line):

    ls -1 *.fastq > list_reads
    ./dsk -file list_reads

## Results visualization

DSK uses a binary format for its output.

By default, the format is HDF5, so you can use HDF5 tools for getting ascii output
from the HDF5 output (such tools are provided with DSK distribution).

In particular, the HDF5 output of DSK holds two sets of data:

* the couples of (kmer,abundance)
* the histogram of abundances

You can get content of a dataset (kmers in partition 0, and whole histogram) with:

    ext/gatb-core/bin/gatb-h5dump -y -d dsk/solid/0          output.h5
    ext/gatb-core/bin/gatb-h5dump -y -d histogram/histogram  output.h5

To see the results as a list of "[kmer] [count]\n", use the `dsk2ascii` binary with the output file from dsk

    dsk2ascii -file output.h5 -out output.txt

To plot kmer coverage distribution,    

    build/ext/gatb-core/bin/gatb-h5dump -y -d histogram/histogram  output.h5  | grep "^\ *[0-9]" | tr -d " " | paste - - | gnuplot -p -e 'plot  "-" with lines'     

### Alternative ways to plot kmer coverage distributions

To plot the kmer coverage distribution, you can also use the `-histo` option and a simple R script provided in the `utils` directory

    ./dsk  -file A.fa -histo 1 -out outputfile
	Rscript utils/plot-histo.R outputfile.histo 
	# the kmer count profile is plotted in the image file outputfile.histo.png
	
With the option `-histo2D`, you can also plot a kmer comparison plot between a genome assembly and a read set to visualize for instance k-mer duplication levels of the given assembly, following a nice idea originally implemented in [KAT (Kmer Analysis Toolkit)](https://github.com/TGAC/KAT) (see examples and how to interpret such plots [in the KAT documentation](http://kat.readthedocs.io/en/latest/walkthrough.html#genome-assembly-analysis-using-k-mer-spectra)). For this purpose, you need to give as input the assembly fasta file and the read files with the option `-file` taking care to give the assembly file as the first file of the list. 

    ./dsk  -file assembly.fa,reads1.fq,reads2.fq -histo2D 1 -out outputname
	# the kmer count matrix is output in file outputname.histo2D, 
	# where lines correspond to abundances in the read dataset and columns abundances in the genome assembly (0 to 10)
	Rscript utils/plot-histo2D.R outputname.histo2D 
	# the kmer comparison plot is plotted in the image file outputname.histo2D.png
	
## Kmers and their reverse complements

DSK converts all kmers to their canonical representation with respect to reverse-complementation.

In other words, a kmer and its reverse complement are considered to be the same object.
For example, with k=3 and assuming the kmer AAA and its reverse complement TTT are both present the input dataset, DSK will consider that one of them is the canonical kmer, e.g. AAA. If AAA is present 2 times and TTT 3 times, then DSK outputs that the count of AAA is 5 and won't return the count of TTT at all.

Note: a canonical kmer is not necessarily the lexicographically smallest one! DSK uses a different ordering for faster performance. Specifically, DSK considers tha A<C<T<G and returns the lexicographically smaller kmer using this alphabet order. So, in the example above, AAA is indeed the canonical kmer.
For the GTA/TAC pair, the lexicographically smallest is GTA however the canonical kmer is TAC (as DSK considers that T<G).


## Larger k-mer sizes

DSK supports arbitrary large k-mer lengths.
Just compile from the source, to support k-mer lengths up to, say, 160, type this in the build folder:

    rm -Rf CMake* && cmake -DKSIZE_LIST="32 64 96 128 160" .. && make

KSIZE_LIST can contain an arbitrary number of multiples of 32.

## Disk space and speed

As a general rule of thumb, run DSK in a folder with plenty of free space, i.e. several times the size of the input dataset.
You can also specify the ```-out-tmp``` parameter to a location with free space.
In the output during execution, the number of passes should be low (below 10). 
Else this generally means that more disk space would make the execution much faster.
DSK auto-detects the free disk space and uses a fraction of it. You can specify ```-max-disk``` parameter to fine-tune this.


## Homepage, contact

https://gatb.inria.fr/software/dsk/

To contact a developer: https://gatb.inria.fr/contact/
