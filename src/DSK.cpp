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

#ifdef OMP
#include <omptl/omptl_algorithm>
#endif

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
const char* DSK::STR_KMER_SIZE          = "-kmer-size";
const char* DSK::STR_MAX_MEMORY         = "-max-memory";
const char* DSK::STR_MAX_DISK           = "-max-disk";
const char* DSK::STR_NKS                = "-nks";
const char* DSK::STR_URI_SOLID_KMERS    = "-solid-kmers";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
DSKAlgorithm<T>::DSKAlgorithm (DSK* dsk)  :
    _bankBinary(0), _kmerSize(0), _nks(0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0),
    _dsk(dsk), _timeInfo(_dsk->_timeInfo)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
DSKAlgorithm<T>::~DSKAlgorithm ()
{
    /** We release the bank handle. */
    if (_bankBinary)  {  delete _bankBinary;  }

    /** We remove physically the partition files. */
    for (size_t i=0; i<_nb_partitions; i++)
    {
        char filename[128];  snprintf (filename, sizeof(filename), getPartitionUri().c_str(), i);
        System::file().remove (filename);
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
template<typename T>
void DSKAlgorithm<T>::execute ()
{
    /** We set the max memory according to the number of used cores. */
    getInput()->setInt (DSK::STR_MAX_MEMORY, getInput()->getInt (DSK::STR_MAX_MEMORY) / getInput()->getInt (DSK::STR_NB_CORES));

    /** Shortcuts attributes. */
    _kmerSize       = getInput()->getInt (DSK::STR_KMER_SIZE);
    _max_memory     = getInput()->getInt (DSK::STR_MAX_MEMORY);
    _max_disk_space = getInput()->getInt (DSK::STR_MAX_DISK);
    _nks            = getInput()->getInt (DSK::STR_NKS);

    /** We add the prefix to the solid kmers uri. */
    getInput()->setStr (DSK::STR_URI_SOLID_KMERS, _dsk->getUriByKey (DSK::STR_URI_SOLID_KMERS) );

    /** We create the binary bank holding the reads in binary format. */
    _bankBinary = new BankBinary (getInput()->getStr (DSK::STR_URI_DATABASE));

    /** We configure dsk by computing the number of passes and partitions we will have
     * according to the allowed disk and memory space. */
    configure ();

    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = _dsk->createIterator<Sequence> (_bankBinary->iterator(), _estimateSeqNb, "DSK");
    LOCAL (itSeq);

    /** We create the solid kmers bag. */
    Bag<T>* solidKmers = createSolidKmersBag ();
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
    getInfo()->add (1, "stats");
    getInfo()->add (2, "solid kmers nb",   "%ld", (System::file().getSize(getInput()->getStr (DSK::STR_URI_SOLID_KMERS)) / sizeof (T)) );
    getInfo()->add (2, "solid kmers uri",  getInput()->getStr (DSK::STR_URI_SOLID_KMERS));

    /** We set the result of the execution. */
    getOutput()->add (0, DSK::STR_URI_SOLID_KMERS,  getInput()->getStr (DSK::STR_URI_SOLID_KMERS));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
void DSKAlgorithm<T>::configure ()
{
    /** We get some information about the bank. */
    _bankBinary->estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));
    u_int64_t bankSize = _estimateSeqTotalSize / MBYTE;

    _volume = kmersNb * sizeof(T) / MBYTE;  // in MBytes

    if (_max_disk_space == 0)  { _max_disk_space = std::min (available_space/2, bankSize);  }
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
    getInfo()->add (1, "config");
    getInfo()->add (2, "current directory", System::file().getCurrentDirectory());
    getInfo()->add (2, "available space",   "%ld", available_space);
    getInfo()->add (2, "bank size",         "%ld", bankSize);
    getInfo()->add (2, "sequence number",   "%ld", _estimateSeqNb);
    getInfo()->add (2, "sequence volume",   "%ld", _estimateSeqTotalSize / MBYTE);
    getInfo()->add (2, "kmers number",      "%ld", kmersNb);
    getInfo()->add (2, "kmers volume",      "%ld", _volume);
    getInfo()->add (2, "max disk space",    "%ld", _max_disk_space);
    getInfo()->add (2, "max memory",        "%ld", _max_memory);
    getInfo()->add (2, "nb passes",         "%d",  _nb_passes);
    getInfo()->add (2, "nb partitions",     "%d",  _nb_partitions);
    getInfo()->add (2, "nb bits per kmer",  "%d",  T::getSize());
    getInfo()->add (2, "nb cores",          "%d",  getDispatcher()->getExecutionUnitsNumber());
}

/********************************************************************************/

template<typename T>
struct Hash  {   T operator () (T& lkmer)
{
    T kmer_hash;

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
template<typename T>
class FillPartitions : public IteratorFunctor
{
public:

    void operator() (Sequence& sequence)
    {
        vector<T> kmers;

        Hash<T> hash;

        /** We build the kmers from the current sequence. */
        model.build (sequence.getData(), kmers);

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            /** We hash the current kmer. */
            T h = hash (kmers[i]);

            /** We check whether this kmer has to be processed during the current pass. */
            if ((h % nbPass) != pass)  { continue; }

            T reduced_kmer = h / nbPass;

            /** We compute in which partition this kmer falls into. */
            size_t p = reduced_kmer % nbPartitions;

            /** We write the kmer into the bag. */
            _cache[p]->insert (kmers[i]);
        }
    }

    FillPartitions (Model<T>& model, size_t nbPasses, size_t currentPass, BagFilePartition<T>& partition)
        : pass(currentPass), nbPass(nbPasses), nbPartitions(partition.size()), _cache(partition,getSynchro()), model(model)
    {
    }

private:
    size_t pass;
    size_t nbPass;
    size_t nbPartitions;
    BagCachePartition<T> _cache;
    Model<T>& model;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
void DSKAlgorithm<T>::fillPartitions (size_t pass, Iterator<Sequence>* itSeq)
{
    TIME_INFO (_timeInfo, "fill partitions");

    /** We create a kmer model. */
    Model<T> model (_kmerSize);

    /** We create the partition files for the current pass. */
    BagFilePartition<T> partitions (_nb_partitions, getPartitionUri());

    /** We launch the iteration of the sequences iterator with the created functors. */
    getDispatcher()->iterate (itSeq, FillPartitions<T> (model, _nb_passes, pass, partitions));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
void DSKAlgorithm<T>::fillSolidKmers (Bag<T>*  solidKmers)
{
    TIME_INFO (_timeInfo, "fill solid kmers");

    Iterator<size_t>* itParts = _dsk->createIterator<size_t> (
        new Range<size_t>::Iterator (0, _nb_partitions-1),
        _nb_partitions,
        "read partitions"
    );
    LOCAL (itParts);

    /** We parse each partition file. */
    for (itParts->first(); !itParts->isDone(); itParts->next())
    {
        /** we build the name of the ith kmers partition file. */
        char filename[128];  snprintf (filename, sizeof(filename), getPartitionUri().c_str(), itParts->item());

        IteratorFile<T> it (filename);
        vector<T> kmers (System::file().getSize(filename) / sizeof(T));

        /** We directly fill the vector from the current partition file. */
        it.fill (kmers);

        /** We check that we got something. */
        if (kmers.empty())  { throw Exception ("DSK: no solid kmers found"); }

#ifdef OMP
        omptl::sort (kmers.begin (), kmers.end ());
#else
        std::sort (kmers.begin (), kmers.end ());
#endif
        u_int32_t max_couv  = 2147483646;
        u_int32_t abundance = 0;
        T previous_kmer = kmers.front();

        for (typename vector<T>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
        {
            T current_kmer = *itKmers;

            if (current_kmer == previous_kmer)  {   abundance++;  }
            else
            {
                /** We should update the abundance histogram => to be done. */

                /** We check that the current abundance is in the correct range. */
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
template<typename T>
Bag<T>* DSKAlgorithm<T>::createSolidKmersBag ()
{
    /** We delete the solid kmers file. */
    System::file().remove (getInput()->getStr (DSK::STR_URI_SOLID_KMERS));

    return new BagCache<T> (new  BagFile<T> (getInput()->getStr (DSK::STR_URI_SOLID_KMERS)), 10*1000);
}

/*********************************************************************************/
/* We define here some concrete subclasses of DSKAlgorithm with specific integer */
/* sizes. According to the end user choice for the kmer size, we can instantiate */
/* one of these classes.                                                         */
/*********************************************************************************/

class DSKAlgorithm32  : public DSKAlgorithm<math::NativeInt64>   {  public: DSKAlgorithm32  (DSK* dsk) : DSKAlgorithm<math::NativeInt64>  (dsk) {}  };
#ifdef INT128_FOUND
class DSKAlgorithm64  : public DSKAlgorithm<math::NativeInt128>  {  public: DSKAlgorithm64  (DSK* dsk) : DSKAlgorithm<math::NativeInt128> (dsk) {}  };
#else
class DSKAlgorithm64  : public DSKAlgorithm<math::LargeInt<2> >  {  public: DSKAlgorithm64  (DSK* dsk) : DSKAlgorithm<math::LargeInt<2> > (dsk) {}  };
#endif
class DSKAlgorithm96  : public DSKAlgorithm<math::LargeInt<3> >  {  public: DSKAlgorithm96  (DSK* dsk) : DSKAlgorithm<math::LargeInt<3> > (dsk) {}  };
class DSKAlgorithm128 : public DSKAlgorithm<math::LargeInt<4> >  {  public: DSKAlgorithm128 (DSK* dsk) : DSKAlgorithm<math::LargeInt<4> > (dsk) {}  };

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
    _parser->add (new OptionOneParam (DSK::STR_KMER_SIZE,       "size of a kmer",                       true            ));
    _parser->add (new OptionOneParam (DSK::STR_MAX_MEMORY,      "max memory in MBytes",                 false,  "1000"  ));
    _parser->add (new OptionOneParam (DSK::STR_MAX_DISK,        "max disk space in MBytes",             false,  "0"     ));
    _parser->add (new OptionOneParam (DSK::STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));
    _parser->add (new OptionOneParam (DSK::STR_URI_SOLID_KMERS, "solid kmers file",                     false,  "solid" ));
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
    size_t kmerSize = _input->getInt (DSK::STR_KMER_SIZE);

    /** According to the kmer size, we instantiate one DSKAlgorithm class and delegate the actual job to it. */
         if (kmerSize < 32)     { DSKAlgorithm32 (this).execute ();  }
    else if (kmerSize < 64)     { DSKAlgorithm64 (this).execute ();  }
    else if (kmerSize < 96)     { DSKAlgorithm96 (this).execute ();  }
    else if (kmerSize < 128)    { DSKAlgorithm128(this).execute ();  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
}
