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

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>

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
const char* DSK::STR_URI_HISTO          = "-histo";

/********************************************************************************/
static const char* progressFormat1 = "DSK: Pass %d/%d, Step 1: partitioning        ";
static const char* progressFormat2 = "DSK: Pass %d/%d, Step 2: counting kmers %2d/%2d";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
DSKAlgorithm<T>::DSKAlgorithm (Tool* dsk)  :
    ToolProxy (dsk),
    _bankBinary(0), _kmerSize(0), _nks(0),
    _progress (0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
    _histogram (0)
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

    setProgress(0);
    if (_histogram)  {  delete _histogram; }
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
    _kmerSize        = getInput()->getInt (DSK::STR_KMER_SIZE);
    _max_memory      = getInput()->getInt (DSK::STR_MAX_MEMORY);
    _max_disk_space  = getInput()->getInt (DSK::STR_MAX_DISK);
    _nks             = getInput()->getInt (DSK::STR_NKS);

    /** We setup the histogram if needed. */
    if (getInput()->get(DSK::STR_URI_HISTO) != 0)  {  _histogram = new Histogram     (10000, getUriByKey (DSK::STR_URI_HISTO));  }
    else                                           {  _histogram = new HistogramNull (); }

    /** We add the prefix to the solid kmers uri. */
    DSK* dsk = dynamic_cast<DSK*> (getRef());   if (!dsk)  { throw Exception ("dynamic cast failed"); }
    getInput()->setStr (DSK::STR_URI_SOLID_KMERS, dsk->getUriByKey (DSK::STR_URI_SOLID_KMERS) );

    /** We create the binary bank holding the reads in binary format. */
    _bankBinary = new BankBinary (getInput()->getStr (DSK::STR_URI_DATABASE));

    /** We configure dsk by computing the number of passes and partitions we will have
     * according to the allowed disk and memory space. */
    configure ();

    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = _bankBinary->iterator();
    LOCAL (itSeq);

    /** We create the solid kmers bag. */
    Bag<T>* solidKmers = createSolidKmersBag ();
    LOCAL (solidKmers);

    /** We configure the progress bar. */
    setProgress ( createIteratorListener (2 * _volume * MBYTE / sizeof(T), "counting kmers"));
    _progress->init ();

    /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
    for (_current_pass=0; _current_pass<_nb_passes; _current_pass++)
    {
        /** 1) We fill the partition files. */
        fillPartitions (_current_pass, itSeq);

        /** 2) We fill the kmers solid file from the partition files. */
        fillSolidKmers (solidKmers);
    }

    _progress->finish ();

    /** We flush the solid kmers file. */
    solidKmers->flush();

    /** We save the histogram if any. */
    _histogram->save ();

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
    if (_volume == 0)   { _volume = 1; }    // tiny files fix

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
        Hash<T>   hash;

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
            _partition[p]->insert (kmers[i]);

            nbWrittenKmers++;
        }

        if (nbWrittenKmers > 50000)   {  _progress.inc (nbWrittenKmers);  nbWrittenKmers = 0;  }
    }

    FillPartitions (Model<T>& model, size_t nbPasses, size_t currentPass, BagFilePartition<T>& partition, IteratorListener* progress)
        : model(model), pass(currentPass), nbPass(nbPasses), nbPartitions(partition.size()), nbWrittenKmers(0),
          _partition(partition,this->newSynchro()), _progress (progress,this->newSynchro())  {}

    ~FillPartitions ()  {}

private:

    Model<T>& model;
    size_t    pass;
    size_t    nbPass;
    size_t    nbPartitions;
    size_t    nbWrittenKmers;

    BagCachePartition<T> _partition;
    ProgressSynchro      _progress;
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
    TIME_INFO (getTimeInfo(), "fill partitions");

    /** We create a kmer model. */
    Model<T> model (_kmerSize);

    /** We create the partition files for the current pass. */
    BagFilePartition<T> partitions (_nb_partitions, getPartitionUri());

    /** We update the message of the progress bar. */
    _progress->setMessage (progressFormat1, _current_pass+1, _nb_passes);

    /** We launch the iteration of the sequences iterator with the created functors. */
    getDispatcher()->iterate (itSeq, FillPartitions<T> (model, _nb_passes, pass, partitions, _progress));
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
    TIME_INFO (getTimeInfo(), "fill solid kmers");

    Iterator<size_t>* itParts = new Range<size_t>::Iterator (0, _nb_partitions-1);
    LOCAL (itParts);

    /** We allocate a vector that will be filled by partitions file (one at time).
     *  We provide an estimation size according to the maximum memory allowed. */
    vector<T> kmers ( (_max_memory * MBYTE) / sizeof(T));

    /** We parse each partition file. */
    for (itParts->first(); !itParts->isDone(); itParts->next())
    {
        /** We update the message of the progress bar. */
        _progress->setMessage (progressFormat2, _current_pass+1, _nb_passes, itParts->item()+1, _nb_partitions);

        /** we build the name of the ith kmers partition file. */
        char filename[128];  snprintf (filename, sizeof(filename), getPartitionUri().c_str(), itParts->item());

        /** We get the length of the current partition file. */
        size_t partitionLen = System::file().getSize(filename) / sizeof(T);

        /** We check that we got something. */
        if (partitionLen == 0)  { throw Exception ("DSK: no solid kmers found"); }

        /** We resize our vector that will be filled with the partition file content.
         * NOTE: we add an extra item and we will set it to the maximum kmer value. */
        kmers.resize (1 + partitionLen);

        /** We directly fill the vector from the current partition file. */
        IteratorFile<T> it (filename);   it.fill (kmers, partitionLen);

        /** We set the extra item to a max value, so we are sure it will sorted at the last location.
         * This trick allows to avoid extra treatment after the loop that computes the kmers abundance. */
        kmers[partitionLen] = ~0;

#ifdef OMP
        omptl::sort (kmers.begin (), kmers.end ());
#else
        std::sort (kmers.begin (), kmers.end ());
#endif
        u_int32_t max_couv  = 2147483646;
        u_int32_t abundance = 0;
        T previous_kmer = kmers.front();

        /** We loop over the sorted solid kmers. */
        for (typename vector<T>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
        {
            if (*itKmers == previous_kmer)  {   abundance++;  }
            else
            {
                /** We should update the abundance histogram*/
                _histogram->inc (abundance);

                /** We check that the current abundance is in the correct range. */
                if (abundance >= _nks  && abundance <= max_couv)
                {
                    solidKmers->insert (previous_kmer);
                }
                abundance     = 1;
                previous_kmer = *itKmers;
            }
        }

        /** We update the progress bar. */
        _progress->inc (kmers.size());
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
    getParser()->add (new OptionOneParam (Tool::STR_URI_DATABASE,   "databank uri",                         true));
    getParser()->add (new OptionOneParam (DSK::STR_KMER_SIZE,       "size of a kmer",                       true            ));
    getParser()->add (new OptionOneParam (DSK::STR_MAX_MEMORY,      "max memory in MBytes",                 false,  "1000"  ));
    getParser()->add (new OptionOneParam (DSK::STR_MAX_DISK,        "max disk space in MBytes",             false,  "0"     ));
    getParser()->add (new OptionOneParam (DSK::STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));
    getParser()->add (new OptionOneParam (DSK::STR_URI_SOLID_KMERS, "solid kmers file",                     false,  "solid" ));
    getParser()->add (new OptionOneParam (DSK::STR_URI_HISTO,       "outputs histogram of kmers abundance", false));
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
    size_t kmerSize = getInput()->getInt (DSK::STR_KMER_SIZE);

    Tool* tool = 0;

    /** According to the kmer size, we instantiate one DSKAlgorithm class and delegate the actual job to it. */
    if (kmerSize < 32)        {  tool = new DSKAlgorithm <math::NativeInt64> (this);  }
    else if (kmerSize < 64)
    {
#ifdef INT128_FOUND
        tool = new DSKAlgorithm <math::NativeInt128> (this);
#else
        tool = new DSKAlgorithm <math::LargeInt<2> > (this);
#endif
    }
    else if (kmerSize < 96)   {  tool = new DSKAlgorithm <math::LargeInt<3> > (this);  }
    else if (kmerSize < 128)  {  tool = new DSKAlgorithm <math::LargeInt<4> > (this);  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }

    /** We use locally the created instance. */
    LOCAL(tool);

    /** We execute the tool. */
    tool->execute ();
}
