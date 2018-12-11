# Copyright 2014 David R. Lamprea.
# Licensed under the European Union Public Licence (EUPL) 1.1.
#
# Exports the following variables:
# * LOOPTOOLS_FOUND
# * LOOPTOOLS_INCLUDE_DIRS - Set of paths to all required headers
# * LOOPTOOLS_LIBRARIES - Set of all required libraries

# The library is called ooptools and not looptools. This is so you can use the -looptools
# option in the compiler.
find_library(LOOPTOOLS_LIBRARIES ooptools HINTS ${LOOPTOOLS}/lib64 ${LOOPTOOLS}/lib ${LOOPTOOLS_LIBDIR})

find_path(LOOPTOOLS_INCLUDE_DIRS clooptools.h HINTS ${LOOPTOOLS}/include ${LOOPTOOLS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LoopTools DEFAULT_MSG LOOPTOOLS_LIBRARIES LOOPTOOLS_INCLUDE_DIRS)
