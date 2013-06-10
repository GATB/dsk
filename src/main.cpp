/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <DSK.hpp>

using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc::impl;
using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
    /** We create an instance of DSK class. */
    DSK dsk;

    /** We create a command line parser for DSK. */
    OptionsParser* parser = dsk.createOptionsParser ();
    LOCAL (parser);    

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We parse the command line arguments. */
        parser->parse (argc, argv);

        /** We execute dsk. */
        IProperties& dskResult = dsk.execute (parser->getProperties());

        /** We may have to dump execution information to stdout. */
        if (dskResult [DSK::STR_QUIET] == 0)
        {
            RawDumpPropertiesVisitor visit;
            dskResult.accept     (&visit);
        }

        /** We may have to dump execution information to stdout. */
        if (dskResult [DSK::STR_STATS_XML] != 0)
        {
            XmlDumpPropertiesVisitor visit (dskResult [DSK::STR_STATS_XML]->getValue());
            dskResult.accept     (&visit);
        }
    }

    catch (OptionFailure& e)
    {
        if (parser->saw("-h"))   {   parser->displayHelp   (stdout);   }
        else                     {   parser->displayErrors (stdout);   }
        return EXIT_FAILURE;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

