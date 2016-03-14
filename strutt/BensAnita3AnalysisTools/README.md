Hello README reader.
If you have doxygen installed then type 

 make docs 

to create a lovely html rendered version of this README at docs/html/index.html.

Alternatively go to http://www.hep.ucl.ac.uk/~strutt/BensAnita3AnalysisTools to see an online version of the documentation.

If neither of those takes your fancy then read on below, ignoring the doxygen tags.

Have fun!

\mainpage 

\tableofcontents

\section intro_sec Introduction

This is the software library I'm using for my analysis.
The idea is to keep any executables doing analysis simple, just a main function that loops over TTrees containing the ANITA-3 data.
The output of each of those executables should be a few summary histograms or a sparsely leaved TTree.
The hard work (turning raw voltage/time data points into interesting numbers) will be done in these libraries.

\section prereq Pre-requisites

\subsection prepre1 1.) ROOT
Available from <a href="https://root.cern.ch/">The ROOT website</a>, if you don't have this, get it pronto.
Compatible with v5.34 and versions greater than 6.

\subsection prepre2 2.) Some of the ANITA offline software libraries (and their pre-requisites)
You will need <a href="http://www.hep.ucl.ac.uk/uhen/anita/eventReader/">EventReaderRoot</a> and <a href="http://www.hep.ucl.ac.uk/uhen/anita/libRootFftwWrapper/">libRootFftwWrapper</a>.

These are, respectively, the fundamental ROOT friendly classes to interact with the ROOT-ified ANITA data, and a set of useful wrapper functions for the brilliant fftw3 library.
Both of these are maintained by <a href="http://www.hep.ucl.ac.uk/~rjn/">Ryan</a>.

\section install_intr Install instructions
Once you have the ANITA environment set up type:
\code{.sh}
svn co https://delos.mps.ohio-state.edu/anitaGround/strutt/BensAnita3AnalysisTools/ BensAnita3AnalysisTools
\endcode
to checkout the latest version from the ANITA SVN.
Then do:
\code{.sh}
make; make install
\endcode
which should install the libraries and header files at the same location as the EventReaderRoot software.

\section probs Foreseeable issues
\subsection makefile_probs Makefile
Compile flags you may need to comment out at various places in the Makefile
\code{.sh}
-wpedantic -std=c++11
\endcode

\section contact Contact
b.strutt.12@ucl.ac.uk