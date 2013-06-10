/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _DSK_HPP_
#define _DSK_HPP_

/********************************************************************************/

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/kmer/impl/Model.hpp>

#include <string>

/********************************************************************************/

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/

/** */
class DSK
{
public:

    /** */
    DSK ();

    /** */
    virtual ~DSK ();
    
    /** */
    IProperties&  execute (IProperties* params);

    /** */
    virtual OptionsParser* createOptionsParser ();

    /** */
    static const char* STR_KMER_SIZE;
    static const char* STR_DB;
    static const char* STR_NB_CORES;
    static const char* STR_MAX_MEMORY;
    static const char* STR_NKS;
    static const char* STR_PREFIX;
    static const char* STR_QUIET;
    static const char* STR_STATS_XML;
    static const char* STR_OUTPUT;

private:

    /** */
    void configure ();
    
    /** */
    void fillPartitions (size_t pass, Iterator<Sequence>* itSeq);

    /** */
    void fillSolidKmers (Bag<kmer_type>*  solidKmers);

    /** */
    virtual Iterator<Sequence>* createSequenceIterator (IteratorListener* progress);

    /** */
    virtual Bag<kmer_type>* createSolidKmersBag ();

    /** */
    void buildBankBinary (Bank& bank);

    /** */
    std::string getPartitionUri ()  {  return _prefix + "%d";  }

    /** */
    std::string getOutputUri ()  { return _prefix + (*_params)[STR_OUTPUT]->getValue(); }

    /** */
    IProperties* _params;
    void setParams (IProperties* params)  { SP_SETATTR(params); }

    /** */
    IProperties* _stats;
    void setStats (IProperties* stats)  { SP_SETATTR(stats); }

    BankBinary* _bankBinary;
    TimeInfo _timeInfo;

    std::string   _filename;
    size_t   _kmerSize;
    size_t   _nks;
    std::string   _prefix;

    ICommandDispatcher* _dispatcher;

    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;
    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;
};

#endif /* _DSK_HPP_ */

