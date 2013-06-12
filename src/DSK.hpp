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

#include <string>

/********************************************************************************/

/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;

/********************************************************************************/

/** \brief Kmer counting class
 */
class DSK : public misc::impl::Tool
{
public:

    /** Constructor. */
    DSK ();

    /** destructor. */
    virtual ~DSK ();
    
    /** */
    static const char* STR_KMER_SIZE;
    static const char* STR_DATABASE;
    static const char* STR_MAX_MEMORY;
    static const char* STR_NKS;
    static const char* STR_PREFIX;
    static const char* STR_SOLID_KMERS;

private:

    /** */
    void  execute ();

    /** */
    void configure ();
    
    /** */
    void fillPartitions (size_t pass, dp::Iterator<Sequence>* itSeq);

    /** */
    void fillSolidKmers (collections::Bag<kmer::impl::kmer_type>*  solidKmers);

    /** */
    virtual collections::Bag<kmer::impl::kmer_type>* createSolidKmersBag ();

    /** */
    void buildBankBinary (IBank& bank);

    /** */
    std::string getPartitionUri ()  {  return _input->getStr (STR_PREFIX) + "%d";  }

    bank::IBank* _bankBinary;

    /** Shortcuts. */
    std::string _filename;
    size_t      _kmerSize;
    size_t      _nks;

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

