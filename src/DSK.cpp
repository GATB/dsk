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
    IBank* bank = Bank::open(props->getStr("-file"));
    
    LOCAL (bank);

    size_t kmerSize = props->getInt (STR_KMER_SIZE);

    size_t abundanceMin = props->getInt (STR_KMER_ABUNDANCE_MIN);
    size_t abundanceMax = props->getInt (STR_KMER_ABUNDANCE_MAX);

    size_t minimizerSize = props->get(STR_MINIMIZER_SIZE)     ? props->getInt(STR_MINIMIZER_SIZE)      : 8;
    size_t minimizerType = props->get(STR_MINIMIZER_TYPE)     ? props->getInt(STR_MINIMIZER_TYPE)      : 0;


    StorageMode_e storageMode = DSK::getStorageMode();

    string output = props->get(STR_URI_OUTPUT) ?
        props->getStr(STR_URI_OUTPUT)   :
        System::file().getBaseName (bank->getId());

    /************************************************************/
    /*                       Storage creation                   */
    /************************************************************/
    Storage* product = StorageFactory(storageMode).create (output, true, false);
    LOCAL (product);

    /************************************************************/
    /*                         Sorting count                    */
    /************************************************************/

    int use_hashing_instead_of_sorting = 0; // 0 = sorting (default)
                                            // 1 = hashing (experimental)

    /** We create a DSK instance and execute it. */
    SortingCountAlgorithm<span> sortingCount (
        product,
        bank, //converter.getResult(),
        kmerSize,
        make_pair(abundanceMin,abundanceMax),
        props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 0,
        props->get(STR_MAX_DISK)   ? props->getInt(STR_MAX_DISK)   : 0,
        props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0,
        gatb::core::tools::misc::KMER_SOLIDITY_DEFAULT,
        props->getInt(STR_HISTOGRAM_MAX),
        use_hashing_instead_of_sorting,
        minimizerType,
        minimizerSize
    );
    sortingCount.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));
    sortingCount.execute();
    dsk.getInfo()->add (1, sortingCount.getInfo());

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
    /** We add options specific to DSK (most important at the end). */
    getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);

    /** We rename the input option. */
    if (IOptionsParser* input = getParser()->getParser (STR_URI_INPUT))  {  input->setName (STR_URI_FILE);  }
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
    else if (kmerSize < KSIZE_4)  { executeAlgorithm <KSIZE_4>  (*this, getInput());  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
}
