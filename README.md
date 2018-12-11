Resummino
=========

**Resummino** computes resummation predictions for Beyond the Standard Model (BSM) particle production at hadron colliders up to the NLO+NLL level. Currently the processes implemented include gaugino-pair production and slepton-pair production. It is able to compute total cross sections as well as invariant-mass and transverse-momentum distributions.

This software is open-source software under the terms of the European Union Public Licence version 1.1 or later.

Currently maintained by David Lamprea and Marcel Rothering at the Research Group of Prof. Dr. Michael Klasen, Institut für Theoretische Physik, Universität Münster, Germany.

Prerequisites
-------------

* Boost (some headers only, required for SLHAea) <http://www.boost.org/>
* GNU Scientific Library (library and headers) <http://www.gnu.org/software/gsl/>
* LHAPDF (library and headers) <http://lhapdf.hepforge.org/>
* QCDLoop (already included; no need to download) <http://qcdloop.fnal.gov/>

Compilation and installation
----------------------------

Download and extract the source tarball and `cd` into it. Then you can use the following commands to compile and install the program:

    $ cmake . [options]
    $ make
    $ make install

Where the possible `[options]` include:

* `-DLHAPDF=/path/to/lhapdf` sets where to find the LHAPDF library, if not in the standard path. The library should be under the `lib` subdirectory (in this case `/path/to/lhapdf/lib`) and the headers should be under `include` (in this case `/path/to/lhapdf/include`). If you want to set the two directories independently you can use `-DLHAPDF_LIB_DIR` and `-DLHAPDF_INCLUDE_DIR` for the library and the headers respectively.
* `-DCMAKE_INSTALL_PREFIX=/path/to/install` sets the path where you want to install **Resummino**.

For further options please consult the **CMake** documentation.

Running the program
-------------------

The process information is stored in an input file. An example of such file can be found in `input/resummino.in`. The SUSY model information should be stored in a separate SUSY Les Houches Accord (SLHA) file (an example `slha.in` is provided), referenced in the input file. Then you can issue:

    $ resummino [options] input_file
    
to run the code. The following options are available:

* `--lo` stops the calculation after the LO result.
* `--nlo` stops the calculation after the NLO result.
* `--parameter-log=params.log` stores the values of all parameters, masses and couplings in a log file `params.log` for future reference.
