#!/bin/bash

# hack to compile DSK faster. if you're not a core GATB dev then you don't need it.
#-rayan

# THIS FILE NEEDS TO BE INTO A SUBFOLDER FROM HERE (e.g. ./build/) via:
# mkdir -p build && cd build && ln -s ../quick-build.sh .

# to compile DSK quickly, do initially:
# cd build && cmake .. && make
# then for each new modification, just do: 
# cd build && ./quick-build.sh

# technical details
#it's massive hack
#it indeed requires that "cmake .. && make" was already run
#also requires the "ccache" software to be installed on the system (http://ccache.samba.org/)
#probably won't be up to date if new files will be added to gatb-core
#debug is untested, but to use it, you also need to have done "cmake -DCMAKE_BUILD_TYPE=Debug .." instead of "cmake .."

mkdir -p obj
rm -f dsk

if [[ "$1" == "-debug" ]]
then
    Debug="Debug"
    debug="_debug"
    debugflag="-g"
    opt="-O0"
else
    if [[ "$1" == "-prof" ]]
    then
        Debug="Debug"
        debug="_debug"
        debugflag="-pg"
        opt="-O0"
    else
        opt="-O3"
    fi
fi

options="$debugflag -I ../../../thirdparty/gatb-core/gatb-core/src/ -I ../src/ \
-I ./ext/gatb-core/include/  \
$opt -DINT128_FOUND  -DWITH_LAMBDA_EXPRESSION  -std=c++0x -Wno-invalid-offsetof "

sources="../src/DSK.cpp ../src/main.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Tool.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/bank/impl/BankBinary.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/bank/impl/BankConverterAlgorithm.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/bank/impl/BankFasta.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/bank/impl/BankRegistery.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/PartitionsCommand.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/SortingCountAlgorithm.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/SortingCountAlgorithmTemplates1.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/SortingCountAlgorithmTemplates2.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/SortingCountAlgorithmTemplates3.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/SortingCountAlgorithmTemplates4.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/Model.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/FileSystemCommon.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/FileSystemLinux.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/FileSystemMacos.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/System.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/SystemInfoCommon.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/ThreadLinux.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/system/impl/ThreadMacos.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/designpattern/impl/Command.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/designpattern/impl/Observer.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Algorithm.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Histogram.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/OptionsParser.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Progress.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/TimeInfo.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Tokenizer.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Property.cpp \
../../../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/XmlReader.cpp 
"

compiler=/usr/local/bin/g++
#compiler=clang++

compile_file(){
    file=$1
    obj=$2
    # so that time command only returns wallclock time on seconds, see http://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
    TIMEFORMAT=%R 
    time=$( { time ccache $compiler -c $options $file -o $obj; } 2>&1 )
    if [[ $time > 1 ]]
    then
        echo "compiling $file took $time secs"
    fi
}

objs=""
for file in $sources 
do
    fname=${file##*/}
    fname=${fname%.*}.o
    objs="$objs obj/$fname"
    compile_file $file obj/$fname &
done

wait

ccache  $compiler $debugflag -o  dsk -rdynamic  -ldl -lpthread -lz -lm $objs ext/gatb-core/lib/$Debug/libhdf5$debug.a
