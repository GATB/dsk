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

#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <omptl/omptl_algorithm>

/********************************************************************************/
// We use the required packages
/********************************************************************************/
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::kmer::impl;

/********************************************************************************/
// We define some string constants.
/********************************************************************************/
const char* DSK::STR_KMER_SIZE   = "-kmer-size";
const char* DSK::STR_DATABASE    = "-db";
const char* DSK::STR_MAX_MEMORY  = "-max-memory";
const char* DSK::STR_NKS         = "-nks";
const char* DSK::STR_PREFIX      = "-prefix";
const char* DSK::STR_SOLID_KMERS = "-solid-kmers";

/********************************************************************************/

struct Hash  {   kmer_type operator () (kmer_type& lkmer)
{
    kmer_type kmer_hash;

    kmer_hash  = lkmer ^ (lkmer >> 14);
    kmer_hash  = (~kmer_hash) + (kmer_hash << 18);
    kmer_hash ^= (kmer_hash >> 31);
    kmer_hash  = kmer_hash * 21;
    kmer_hash ^= (kmer_hash >> 11);
    kmer_hash += (kmer_hash << 6);
    kmer_hash ^= (kmer_hash >> 22);

    return kmer_hash;
}};

/********************************************************************************/
class FillPartitions
{
public:

    void operator() (Sequence& sequence)
    {
        vector<kmer_type> kmers;

        Hash hash;

        /** We build the kmers from the current sequence. */
        model.build (sequence.getData(), kmers);

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            /** We hash the current kmer. */
            kmer_type h = hash (kmers[i]);

            /** We check whether this kmer has to be processed during the current pass. */
            if ((h % nbPass) != pass)  { continue; }

            kmer_type reduced_kmer = h / nbPass;

            /** We compute in which partition this kmer falls into. */
            size_t p = reduced_kmer % _cache.size();

            /** We write the kmer into the bag. */
            _cache[p]->insert (kmers[i]);
        }
    }

    FillPartitions (KmerModel& model, size_t nbPasses, size_t currentPass, BagFilePartition<kmer_type>& partition, ISynchronizer* synchro)
        : pass(currentPass), nbPass(nbPasses), _cache(partition,synchro), model(model)
    {
    }

private:
    size_t pass;
    size_t nbPass;
    BagCachePartition<kmer_type> _cache;
    KmerModel& model;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
DSK::DSK ()
    : Tool ("dsk"), _bankBinary(0), _filename(""), _kmerSize(27), _nks(3),
      _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
      _max_disk_space(0), _max_memory(1000), _volume(0), _nb_passes(0), _nb_partitions(0)
{
    /** We add options specific to DSK. */
    _parser->add (new OptionOneParam (DSK::STR_KMER_SIZE,   "size of a kmer",                       true            ));
    _parser->add (new OptionOneParam (DSK::STR_DATABASE,    "URI of the bank",                      true            ));
    _parser->add (new OptionOneParam (DSK::STR_MAX_MEMORY,  "max memory",                           false,  "1000"  ));
    _parser->add (new OptionOneParam (DSK::STR_NKS,         "abundance threshold for solid kmers",  false,  "3"     ));
    _parser->add (new OptionOneParam (DSK::STR_SOLID_KMERS, "solid kmers file",                     false,  "solid" ));
    _parser->add (new OptionOneParam (DSK::STR_PREFIX,      "prefix URI for temporary files",       false,  "dsk."  ));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
DSK::~DSK ()
{
    if (_bankBinary)  {  delete _bankBinary;  }
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
    _kmerSize   = _input->getInt (STR_KMER_SIZE);
    _filename   = _input->getStr (STR_DATABASE);
    _max_memory = _input->getInt (STR_MAX_MEMORY);
    _nks        = _input->getInt (STR_NKS);

    /** We add the prefix to the solid kmers uri. */
    _input->setStr (STR_SOLID_KMERS, _input->getStr (STR_PREFIX) + _input->getStr (STR_SOLID_KMERS));

    /** We create the binary bank holding the reads in binary format. */
    _bankBinary = new BankBinary (_filename + ".bin");

    /** We configure dsk. */
    configure ();

    // We create the sequences iterator.
    Iterator<Sequence>* itSeq = createIterator<Sequence> (_bankBinary->iterator(), _estimateSeqNb, "DSK");
    LOCAL (itSeq);

    // We create the solid kmers bag
    Bag<kmer_type>* solidKmers = createSolidKmersBag ();
    LOCAL (solidKmers);

    /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
    for (size_t pass=0; pass<_nb_passes; pass++)
    {
        /** 1) We fill the partition files. */
        fillPartitions (pass, itSeq);

        /** 2) We fill the kmers solid file from the partition files. */
        fillSolidKmers (solidKmers);
    }

    /** We flush the solid kmers file. */
    solidKmers->flush();

    /** We gather some statistics. */
    _info->add (1, "stats");
    _info->add (2, "solid kmers nb",   "%ld", (System::file().getSize(_input->getStr (STR_SOLID_KMERS)) / sizeof (kmer_type)) );
    _info->add (2, "solid kmers uri",  _input->getStr (STR_SOLID_KMERS));

    /** We set the result of the execution. */
    _output->add (0, STR_SOLID_KMERS,  _input->getStr (STR_SOLID_KMERS));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void DSK::configure ()
{
    // We create a Bank instance.
    Bank bank (_filename);

    // We get some estimations about the bank
    bank.estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    u_int64_t bankSize = bank.getSize() / MBYTE;
    u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));

    _volume = kmersNb * sizeof(kmer_type) / MBYTE;  // in MBytes

    _max_disk_space = std::min (available_space/2, bankSize);

    if (_max_disk_space == 0)  { _max_disk_space = 10000; }

    _nb_passes = ( _volume / _max_disk_space ) + 1;

    size_t max_open_files = System::file().getMaxFilesNumber() / 2;
    u_int64_t volume_per_pass;

    do  {
        volume_per_pass = _volume / _nb_passes;
        _nb_partitions  = ( volume_per_pass / _max_memory ) + 1;

        if (_nb_partitions >= max_open_files)   { _nb_passes++;  }
        else                                    { break;         }

    } while (1);

    /** We gather some statistics. */
    _info->add (1, "config");
    _info->add (2, "current directory", System::file().getCurrentDirectory());
    _info->add (2, "available space",   "%ld", available_space);
    _info->add (2, "bank size",         "%ld", bankSize);
    _info->add (2, "sequence number",   "%ld", _estimateSeqNb);
    _info->add (2, "sequence volume",   "%ld", _estimateSeqTotalSize / MBYTE);
    _info->add (2, "kmers number",      "%ld", kmersNb);
    _info->add (2, "kmers volume",      "%ld", _volume);
    _info->add (2, "max disk space",    "%ld", _max_disk_space);
    _info->add (2, "nb passes",         "%d",  _nb_passes);
    _info->add (2, "nb partitions",     "%d",  _nb_partitions);
    _info->add (2, "nb bits per kmer",  "%d",  kmer_type::getSize());
    _info->add (2, "nb cores",          "%d",  _dispatcher->getExecutionUnitsNumber());

    // We may have to build the binary bank if not already existing.
    buildBankBinary (bank);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void DSK::fillPartitions (size_t pass, Iterator<Sequence>* itSeq)
{
    TIME_INFO (_timeInfo, "fill partitions");

    /** We create a kmer model. */
    KmerModel model (_kmerSize);

    /** We create the partition files for the current pass. */
    BagFilePartition<kmer_type> partitions (_nb_partitions, getPartitionUri().c_str());

    /** We create a shared synchronizer for the partitions building. */
    ISynchronizer* synchro = System::thread().newSynchronizer();

    /** We launch the iteration of the sequences iterator with the created functors. */
    _dispatcher->iterate (*itSeq, FillPartitions (model, _nb_passes, pass, partitions, synchro));

    /** We cleanup resources. */
    delete synchro;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void DSK::fillSolidKmers (Bag<kmer_type>*  solidKmers)
{
    TIME_INFO (_timeInfo, "fill solid kmers");

    /** We parse each partition file. */
    for (size_t i=0; i<_nb_partitions; i++)
    {
        char filename[128];  snprintf (filename, sizeof(filename), getPartitionUri().c_str(), i);

        IteratorFile<kmer_type> it (filename);
        vector<kmer_type> kmers;

        it.fill (kmers);

#ifdef OMP
        omptl::sort (kmers.begin (), kmers.end ());
#else
        std::sort (kmers.begin (), kmers.end ());
#endif

        u_int32_t max_couv  = 2147483646;
        u_int32_t abundance = 0;
        kmer_type previous_kmer = kmers.front();

        for (vector<kmer_type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
        {
            kmer_type current_kmer = *itKmers;

            if (current_kmer == previous_kmer)  {   abundance++;  }
            else
            {
                if (abundance >= _nks  && abundance <= max_couv)
                {
                    solidKmers->insert (previous_kmer);
                }
                abundance = 1;
            }
            previous_kmer = current_kmer;
        }

        if (abundance >= _nks && abundance <= max_couv)
        {
            solidKmers->insert (previous_kmer);
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bag<kmer_type>* DSK::createSolidKmersBag ()
{
    /** We delete the solid kmers file. */
    System::file().remove (_input->getStr (STR_SOLID_KMERS));

    return new BagCache<kmer_type> (new  BagFile<kmer_type> (_input->getStr (STR_SOLID_KMERS)), 5*1000);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void DSK::buildBankBinary (IBank& bank)
{
    TIME_INFO (_timeInfo, "bank conversion");

    Progress* progress = 0;

    if ( _input->get(STR_QUIET) == 0)  {  progress =  new Progress (_estimateSeqNb, "FASTA to binary conversion");  }

    // We convert the FASTA bank in binary format
    IProperties* props = BankHelper::singleton().convert (bank, *_bankBinary, progress);
    LOCAL (props);

    _info->add (1, props);
}
 
