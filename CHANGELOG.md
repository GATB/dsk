# Change Log

--------------------------------------------------------------------------------
## [Unreleased]

--------------------------------------------------------------------------------
## [2.3.0] - 2018-06-13

Many new features and bug fixes since latest official release, including:

* Switching gatb-core version from 1.2.1 to 1.4.1 (see the details in file `thirdparty/gatb-core/gatb-core/RELEASES.md` or on the [github releases page](https://github.com/GATB/gatb-core/releases)). Importantly for DSK :
	
	* Faster k-mer counting (inspired by KMC3 but not yet as fast :).
	* Compiling GATB-Core library now requires c++/11 capable compilers.
	* CMake 3.1.0 is the minimum release of CMake required to compile GATB-Core.
	* Bug fixes in some multi-threaded situations.

* Some new features (also coming from updates in gatb-core but not yet in an official gatb-core release):

	* An easier way to plot the kmer abundance profile (option `-histo` and R script in the `utils/` directory).
	* New feature: can output kmer count abundance matrix between a genome assembly and sequencing read datasets to plot a kmer comparison plot, inspired by [KAT (Kmer Analysis Toolkit)](https://github.com/TGAC/KAT) (option `-histo2D` and R script in the `utils/` directory).
	* New custom solidity option to output kmers specific to a subset of the input files. (`-solidity-custom`).


--------------------------------------------------------------------------------
## [2.2.2] - 2016-07-16

This is a bugfix release since 2.1.0. No new features. More stability.

Note: in the binary release tar file, the correct `simple_test.sh` script is in the `test/` folder.

Thanks to Hamid Mohamadi for noting a mistake in the README.

--------------------------------------------------------------------------------
## [1.0.6] - 2016-03-25

Initial github release. disregard the version name, this is DSK v2 not v1.