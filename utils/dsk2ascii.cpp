//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <iostream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                       Dump solid kmers in ASCII format                       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides a graph URL (supposed to be in HDF5 format).
    if (argc < 2)
    {
        std::cerr << "You must provide a HDF5 file." << std::endl;
        return EXIT_FAILURE;
    }

    // We get a handle on the HDF5 storage object.
    // Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
    Storage* storage = StorageFactory(STORAGE_HDF5).load (argv[1]);
    LOCAL (storage);

    // We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
    Collection<Kmer<>::Count>& solidKmers = storage->getGroup("dsk").getCollection<Kmer<>::Count> ("solid");

    // We can access each of these information through a Properties object
    Properties props;   props.readXML (solidKmers.getProperty("properties"));

    // We create a Model instance. It will help to dump the kmers in
    // a human readable form (ie as a string of nucleotides)
    Kmer<>::Model model (props.getInt ("kmer_size"));

    // We iterate (through a lambda expression) the solid kmers from the retrieved collection
    solidKmers.iterate ([&] (const Kmer<>::Count& count)
    {
        cout << model.toString(count.value) << " " << count.abundance << endl;
    });
}
//! [snippet1]
