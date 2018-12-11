# Copyright 2014 David R. Lamprea.
# Licensed under the European Union Public Licence (EUPL) 1.1.
#
# Exports the following variables:
# * GSL_FOUND
# * GSL_INCLUDE_DIRS - Set of paths to all required headers
# * GSL_LIBRARIES - Set of all required libraries

# The required headers are:
#  gsl/gsl_blas.h
#  gsl/gsl_math.h
#  gsl/gsl_vector.h
#  gsl/gsl_matrix.h
#  gsl/gsl_integration.h
#  gsl/gsl_monte_vegas.h
#  gsl/gsl_multifit_nlin.h

find_path(GSL_INCLUDE_DIR gsl/gsl_math.h)

find_library(GSL_LIBRARY gsl)
find_library(GSL_CBLAS_LIBRARY NAMES cblas gslcblas)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARY GSL_CBLAS_LIBRARY GSL_INCLUDE_DIR)

if (GSL_FOUND)
  set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
  set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
endif ()

mark_as_advanced(GSL_LIBRARY GSL_CBLAS_LIBRARY GSL_INCLUDE_DIR)
