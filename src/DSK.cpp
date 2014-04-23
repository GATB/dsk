/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
    IBank* bank = BankRegistery::singleton().getFactory()->createBank(props->getStr("-file"));
    
    LOCAL (bank);

    size_t kmerSize = props->get(STR_KMER_SIZE)  ? props->getInt(STR_KMER_SIZE) : 31;
    size_t nks      = props->get(STR_NKS)        ? props->getInt(STR_NKS)       : 3;

    StorageMode_e storageMode = DSK::getStorageMode();

    string output = props->get(STR_URI_OUTPUT) ?
        props->getStr(STR_URI_OUTPUT)   :
        System::file().getBaseName (bank->getId());

    string binaryBankUri = System::file().getCurrentDirectory() + "/bank.bin";

    /************************************************************/
    /*                       Storage creation                   */
    /************************************************************/
    Storage* product = StorageFactory(storageMode).create (output, true, false);
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
        props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 0,
        props->get(STR_MAX_DISK)   ? props->getInt(STR_MAX_DISK)   : 0,
        props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
    );
    sortingCount.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));
    sortingCount.execute();
    dsk.getInfo()->add (1, sortingCount.getInfo());

    /** We can get rid of the binary bank. */
    System::file().remove (binaryBankUri);

    /** We set the output uri. */
    dsk.getOutput()->add (0, STR_KMER_SOLID,  output);
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
    /** We get an OptionsParser for DSK. */
    OptionsParser parser = getOptionsParser();

    /** We add options specific to DSK (most important at the end). */
    getParser()->push_front (parser);
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
         if (kmerSize < KSIZE_1)  { executeAlgorithm <KSIZE_1>  (*this, getInput());  }
    else if (kmerSize < KSIZE_2)  { executeAlgorithm <KSIZE_2>  (*this, getInput());  }
    else if (kmerSize < KSIZE_3)  { executeAlgorithm <KSIZE_3>  (*this, getInput());  }
    else if (kmerSize < KSIZE_4)  { executeAlgorithm <KSIZE_4> (*this, getInput());  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
OptionsParser DSK::getOptionsParser (bool includeMandatory)
{
    OptionsParser parser ("DSK");

    /** We add options specific to DSK (most important at the end). */
    parser.push_front (new OptionOneParam (STR_URI_FILE,        "file containing reads (e.g. FASTA/FASTQ)",true));
    parser.push_front (new OptionOneParam (STR_KMER_SIZE,       "size of a kmer",                       true            ));
    parser.push_front (new OptionOneParam (STR_URI_OUTPUT,      "output file (if not set basename of the input file)", false));
    parser.push_front (new OptionOneParam (STR_MAX_MEMORY,      "max memory in MBytes",                 false,  "1000"  ));
    parser.push_front (new OptionOneParam (STR_MAX_DISK,        "max disk space in MBytes",             false,  "0"     ));
    parser.push_front (new OptionOneParam (STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));

    return parser;
}
