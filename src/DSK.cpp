/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <DSK.hpp>

using namespace std;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
static void executeAlgorithm (DSK& dsk, IProperties* props)
{
    IBank* bank = new BankFasta (props->getStr(STR_URI_DB));
    LOCAL (bank);

    size_t kmerSize = props->get(STR_KMER_SIZE)  ? props->getInt(STR_KMER_SIZE) : 31;
    size_t nks      = props->get(STR_NKS)        ? props->getInt(STR_NKS)       : 3;

    StorageMode_e storageMode =  props->getInt(STR_OUTPUT_FORMAT) == 0 ?  STORAGE_FILE : STORAGE_HDF5;

    string output = props->get(STR_URI_OUTPUT) ?
        props->getStr(STR_URI_OUTPUT)   :
        System::file().getBaseName (bank->getId());

    string binaryBankUri = System::file().getCurrentDirectory() + "/bank.bin";

    /************************************************************/
    /*                       Storage creation                   */
    /************************************************************/
    Storage* product = StorageFactory(storageMode).createStorage (output, true, false);
    LOCAL (product);

    /************************************************************/
    /*                         Bank conversion                  */
    /************************************************************/
    /** We create the binary bank. */
    BankConverterAlgorithm converter (bank, kmerSize, binaryBankUri);
    converter.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));
    converter.execute();
    dsk.getInfo()->add (1, converter.getInfo());

    /************************************************************/
    /*                         Sorting count                    */
    /************************************************************/
    /** We create a DSK instance and execute it. */
    SortingCountAlgorithm<span> sortingCount (
        product,
        converter.getResult(),
        kmerSize,
        nks,
        props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 1000,
        props->get(STR_MAX_DISK)   ? props->getInt(STR_MAX_DISK)   : 0,
        props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
    );
    sortingCount.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));
    sortingCount.execute();
    dsk.getInfo()->add (1, sortingCount.getInfo());

    /** We can get rid of the binary bank. */
    System::file().remove (binaryBankUri);

    /** We set the output uri. */
    dsk.getOutput()->add (0, STR_KMER_SOLID,  props->getStr(STR_PREFIX) + dsk.getInput()->getStr (STR_KMER_SOLID));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
DSK::DSK () : Tool ("dsk")
{
    /** We add options specific to DSK (most important at the end). */
    getParser()->push_front (new OptionOneParam (STR_VERBOSE,         "verbosity level",                      false,  "1"));
    getParser()->push_front (new OptionOneParam (STR_OUTPUT_FORMAT,   "output format (0 binary, 1 HDF5)",     false,  "0"));
    getParser()->push_front (new OptionOneParam (STR_PREFIX,          "prefix for output files",              false, "tmp."));
    getParser()->push_front (new OptionOneParam (STR_PARTITION_TYPE,  "partitioning type : 0 for map (default), 1 for vector", false, "0"));
    getParser()->push_front (new OptionOneParam (STR_URI_HISTOGRAM,   "outputs histogram of kmers abundance", false));
    getParser()->push_front (new OptionOneParam (STR_KMER_SOLID,      "solid kmers file",                     false,  "solid" ));
    getParser()->push_front (new OptionOneParam (STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));
    getParser()->push_front (new OptionOneParam (STR_MAX_DISK,        "max disk space in MBytes",             false,  "0"     ));
    getParser()->push_front (new OptionOneParam (STR_MAX_MEMORY,      "max memory in MBytes",                 false,  "1000"  ));
    getParser()->push_front (new OptionOneParam (STR_KMER_SIZE,       "size of a kmer",                       true            ));
    getParser()->push_front (new OptionOneParam (STR_URI_DB,          "databank uri",                         true));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void DSK::execute ()
{
    /** we get the kmer size chosen by the end user. */
    size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    /** According to the kmer size, we instantiate one DSKAlgorithm class and delegate the actual job to it. */
         if (kmerSize < 32)   { executeAlgorithm <32>  (*this, getInput());  }
    else if (kmerSize < 64)   { executeAlgorithm <64>  (*this, getInput());  }
    else if (kmerSize < 96)   { executeAlgorithm <96>  (*this, getInput());  }
    else if (kmerSize < 128)  { executeAlgorithm <128> (*this, getInput());  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
}
