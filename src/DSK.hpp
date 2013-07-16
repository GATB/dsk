/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _DSK_HPP_
#define _DSK_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>

#include <string>

/********************************************************************************/

/** NOTE: we should not include namespaces here => only to make developer life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;

/********************************************************************************/

/** \brief Kmer counting class
 *
 * This is the high level class for running DSK, the counting kmer algorithm. It is an
 * implementation of the Tool interface.
 *
 * The real job is delegated to another class: DSKAlgorithm.
 * Actually, DSK analyzes the kmer size chosen by the user and
 * then decides what kind of DSKAlgorithm to use according to the
 * kmer size.
 */
class DSK : public misc::impl::Tool
{
public:

    /** Constructor. */
    DSK ();

    /** */
    static const char* STR_KMER_SIZE;
    static const char* STR_MAX_MEMORY;
    static const char* STR_MAX_DISK;
    static const char* STR_NKS;
    static const char* STR_URI_SOLID_KMERS;
    static const char* STR_URI_HISTO;
    static const char* STR_PARTITION_TYPE;

private:

    /** \copydoc Tool::execute. */
    void  execute ();

    /** The "real" class DSKAlgorithm needs to know the input parameters only
     * provided by DSK class (as a subclass of Tool). */
    template<typename T> friend class DSKAlgorithm;
};

/********************************************************************************/

/** \brief Class performing the DSK counting.
 *
 * This class does the real job of counting the kmers from a reads database.
 *
 * Note that it is intended to be used as a delegate by the DSK class.
 *
 * This is a template class whose template argument is the kind of integer used for
 * kmers (integers on 64 bits, 128 bits, etc...)
 *
 * We define some template instantiations of this DSKAlgorithm; such an instantiation
 * does the real job of kmers counting. By defining several instantiations, we allow
 * to choose dynamically the correct class according to the user choice for kmer size
 * (remember that initial Minia version had to be re-compiled for different kmer size).
 */
template<typename T> class DSKAlgorithm : public misc::impl::ToolProxy
{
public:

    /** Constructor.
     * \param[in] dsk : Tool instance for dsk that holds input parameters. */
    DSKAlgorithm (misc::impl::Tool* dsk);

    /** Destructor */
    virtual ~DSKAlgorithm ();

    /** Compute several values, in particular the number of passes and partitions. */
    void configure (bank::IBank* bank);

    /** Fill partition files (for a given pass) from a sequence iterator.
     * \param[in] pass  : current pass whose value is used for choosing the partition file
     * \param[in] itSeq : sequences iterator whose sequence are cut into kmers to be split.
     */
    void fillPartitions (size_t pass, dp::Iterator<Sequence>* itSeq);

    /** Fill the solid kmers bag from the partition files (one partition after another one).
     * \param[in] solidKmers : bag to put the solid kmers into.
     */
    void fillSolidKmers (collections::Bag<T>*  solidKmers);

    /** */
    std::vector <size_t> getNbCoresList ();

    /** Create the bag in which the solid kmers will be put into. The actual kind of bag is likely
     * to be a file.
     * \return the solid kmers bag.
     */
    virtual collections::Bag<T>* createSolidKmersBag ();

    /** Compute the format of the URI for the partition files. This format is set from the user preferences,
     * in particular a suffix may be used.
     * \return the output format.
     */
    std::string getPartitionFormat () {  return getInput()->getStr (DSK::STR_URI_PREFIX) + "partition.%d";  }


    std::string getPartitionUri (size_t i)
    {
        char filename[128];  snprintf (filename, sizeof(filename), getPartitionFormat().c_str(), i);
        return filename;
    }

    /** Process the kmers counting. It is mainly composed of a loop over the passes, and for each pass
     * 1) we build the partition files then 2) we fill the solid kmers file from the partitions.
     */
    void  execute ();

    /** Get the iterator listener instance. */
    dp::IteratorListener* getProgress ()  { return _progress; }

    /** */
    size_t getNks ()  { return _nks; }

private:

    /** Shortcuts for the user input parameters. . */
    size_t      _kmerSize;
    size_t      _nks;
    size_t      _partitionType;
    size_t      _nbCores;

    dp::IteratorListener* _progress;
    void setProgress (dp::IteratorListener* progress)  { SP_SETATTR(progress); }

    /** Values computed for algorithm parameterization. In particular, we have one value for the number
     * of passes and one value for the number of partitions.
     * Such values are computed both:
     *      - from system resources (file system resources, memory resources)
     *      - user preferences (max disk space, max memory)
     */
    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;
    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;
    u_int32_t _current_pass;

    misc::impl::Histogram* _histogram;
};

/********************************************************************************/

#endif /* _DSK_HPP_ */

