
--------------------------------------------------------------------------------
How to provide several files to DSK ?
--------------------------------------------------------------------------------

If your bank is split in several files A1,A2,...,An, you can provide all the files to DSK by 
providing the list of files to the '-file' option. EX:
    dsk  -file A1.fa,A2.fa,A3.fa  -kmer-size 31
    

--------------------------------------------------------------------------------
How to get textual output from the DSK output ?
--------------------------------------------------------------------------------

DSK uses a binary format for its output. 

By default, the format is HDF5, so you can use HDF5 tools for getting ascii output
from the HDF5 output (such tools are provided with DSK distribution).

In particular, the HDF5 output of DSK holds two sets of data:
    - the couples of (kmer,abundance)
    - the histogram of abundances
    
You can get a quick look at the datasets with:
    h5dump -n output.h5

If you want to dump the histogram content in textual form, you can do the following:
    h5dump -y -d dsk/histogram  output.h5
-> you will get all the information of the histogram dataset.

For instance, you can get a gnuplot output of the histogram dataset:
    h5dump -y -d dsk/histogram  output.h5 | grep [0-9] | grep -v [A-Z].* | paste - - | gnuplot -p -e 'plot  "-" with lines'     

If you are interested in getting only the ascii values of a dataset in a 'values.txt' file, you can do:
    h5dump -y -d dsk/solid -o values.txt output.h5


