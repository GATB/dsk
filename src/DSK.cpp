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

/********************************************************************************/

#if 1
#define ProductFactoryCurrent ProductFileFactory
#else
#define ProductFactoryCurrent ProductHDF5Factory
#endif

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
static void executeAlgorithm (DSK& dsk, IProperties* props)
{
    IBank* bank = new Bank (props->getStr(STR_URI_DB));
    LOCAL (bank);

    size_t kmerSize = props->get(STR_KMER_SIZE)  ? props->getInt(STR_KMER_SIZE) : 27;
    size_t nks      = props->get(STR_NKS)        ? props->getInt(STR_NKS)       : 3;
    string bargraph = props->get(STR_VERBOSE) ?  "2" : "0";

    string output = props->get(STR_URI_OUTPUT) ?
        props->getStr(STR_URI_OUTPUT)   :
        System::file().getBaseName (bank->getId());

    string binaryBankUri = System::file().getCurrentDirectory() + "/bank.bin";

    /************************************************************/
    /*                       Product creation                   */
    /************************************************************/
    Product<ProductFactoryCurrent>* product = ProductFactoryCurrent::createProduct (output, true, false);
    LOCAL (product);

    /************************************************************/
    /*                         Bank conversion                  */
    /************************************************************/
    /** We create the binary bank. */
    BankConverterAlgorithm converter (bank, kmerSize, binaryBankUri);
    converter.getInput()->add (0, STR_PROGRESS_BAR, bargraph);
    converter.execute();
    dsk.getInfo()->add (1, converter.getInfo());

    /************************************************************/
    /*                         Sorting count                    */
    /************************************************************/
    /** We create a DSK instance and execute it. */
    SortingCountAlgorithm<ProductFactoryCurrent, T> sortingCount (
        product,
        converter.getResult(),
        kmerSize,
        nks,
        props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 1000,
        props->get(STR_MAX_DISK)   ? props->getInt(STR_MAX_DISK)   : 0,
        props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
    );
    sortingCount.getInput()->add (0, STR_PROGRESS_BAR, bargraph);
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
    /** We add options specific to DSK. */
    getParser()->add (new OptionOneParam (STR_URI_DB,          "databank uri",                         true));
    getParser()->add (new OptionOneParam (STR_KMER_SIZE,       "size of a kmer",                       true            ));
    getParser()->add (new OptionOneParam (STR_MAX_MEMORY,      "max memory in MBytes",                 false,  "1000"  ));
    getParser()->add (new OptionOneParam (STR_MAX_DISK,        "max disk space in MBytes",             false,  "0"     ));
    getParser()->add (new OptionOneParam (STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));
    getParser()->add (new OptionOneParam (STR_KMER_SOLID,      "solid kmers file",                     false,  "solid" ));
    getParser()->add (new OptionOneParam (STR_URI_HISTOGRAM,   "outputs histogram of kmers abundance", false));
    getParser()->add (new OptionOneParam (STR_PARTITION_TYPE,  "partitioning type : 0 for map (default), 1 for vector", false, "0"));
    getParser()->add (new OptionOneParam (STR_PREFIX,          "prefix for output files",              false, "tmp."));
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
         if (kmerSize < 32)   { executeAlgorithm <LargeInt<1> > (*this, getInput());  }
    else if (kmerSize < 64)   { executeAlgorithm <LargeInt<2> > (*this, getInput());  }
    else if (kmerSize < 96)   { executeAlgorithm <LargeInt<3> > (*this, getInput());  }
    else if (kmerSize < 128)  { executeAlgorithm <LargeInt<4> > (*this, getInput());  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
}
