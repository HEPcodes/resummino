Resummino
=========

Resummino computes predictions for selected Beyond the Standard Model (BSM) processes at hadron colliders.
In the framework of the Minimal Supersymmetric Standard Model (MSSM) predictions at NNLO+NNLL are available for slepton-pair and electroweakino-pair production (these processes are also available with a collinear improved version of threshold resummation for NLO+NLL); predictions at LO are available for the associated gaugino-gluino and gaugino-squark production.
In general BSM models with extended gauge sectors, predictions for Z' and W' production are available at approximate NNLO+NNLL accuracy.
The code is able to compute total cross sections as well as invariant-mass and transverse-momentum distributions (only for Drell-Yan-like processes). 
This software is open-source software under the terms of the European Union Public Licence version 1.1 or later.

Currently maintained by Juri Fiaschi in the Particle Physics Group at the University of Liverpool and Alexander Neuwirth in the Research Group of Prof. Dr. Michael Klasen at the Institut für Theoretische Physik, Universität Münster, Germany.

Prerequisites
-------------

* Boost (some headers only; required for SLHAea) <http://www.boost.org/>
* GNU Scientific Library <http://www.gnu.org/software/gsl/>
* LHAPDF <http://lhapdf.hepforge.org/>
* LoopTools <http://www.feynarts.de/looptools//>

Compilation and installation
----------------------------

Download and extract the source tarball and `cd` into it. Then you can use the following commands to compile and install the program:

    $ cmake . [options]
    $ make
    $ make install

Where the possible `[options]` include:

* `-DLHAPDF=/path/to/lhapdf` sets where to find the LHAPDF library, if not in the standard path. The library should be under the `lib` subdirectory (in this case `/path/to/lhapdf/lib`) and the headers should be under `include` (in this case `/path/to/lhapdf/include`). If you want to set the two directories independently you can use `-DLHAPDF_LIB_DIR` and `-DLHAPDF_INCLUDE_DIR` for the library and the headers respectively.
* `-DLOOPTOOLS=/path/to/looptools` sets where to find the LoopTools library, if not in the standard path. The library should be in a `lib` or `lib64` subdirectory and the headers should be under `include`. For LoopTools-2.13 you usually have to use for instance `-DLOOPTOOLS=~/.lib/LoopTools-2.13/x86_64-Linux` if you have installed the library in a folder ~/.lib.
* `-DCMAKE_PREFIX_PATH=path/to/gsllib` sets the path where you have installed gsl, if not in the standard path.
* `-DCMAKE_INSTALL_PREFIX=/path/to/install` sets the path where you want to install Resummino.

For further options please consult the CMake documentation.

Running the program
-------------------

The process information is stored in an input file. An example of such file can be found in `input/resummino.in`. This input files should reference two model files: The SUSY model information should be stored in a separate SUSY Les Houches Accord (SLHA) file (an example `slha.in` is provided), and the new gauge bosons model should be stored in an input file with a similar format (an example defining the Sequential Standard Model is stored in `input/ssm.in`).
NOTE: an `slha.in` file should always be present, even when computing Z'/W' cross sections and SUSY particle masses would not be considered.

To run the code:

    $ resummino [options] input_file
    
The following options are available:

* `--lo` stops the calculation after the LO result.
* `--nlo` stops the calculation after the NLO result.
* `--nll` computes ordinary threshold resummation at NLO+NLL.
* `[no flag]` computes ordinary threshold resummation at NLO+NLL; when available collinear improved resummation at NLO+NLL is used.
* `--nnll` computes ordinary threshold resummation at the approximate NNLO+NNLL accuracy level for Drell-Yan like processes (available for lepton-pair, slepton-pair and elecroweakino-pair production).
* `--parameter-log=params.log` stores the values of all parameters, masses and couplings in a log file `params.log` for future reference.
