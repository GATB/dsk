/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _DSK_HPP_
#define _DSK_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

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
class DSK : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    DSK ();

private:

    /** \copydoc Tool::execute. */
    void  execute ();
};

/********************************************************************************/

#endif /* _DSK_HPP_ */

