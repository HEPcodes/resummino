set(BOOST_REQUIRED_HEADERS
  boost/algorithm/string/classification.hpp
  boost/algorithm/string/join.hpp
  boost/algorithm/string/predicate.hpp
  boost/algorithm/string/split.hpp
  boost/lexical_cast.hpp
)

foreach(headername ${BOOST_REQUIRED_HEADERS})
  string(TOUPPER "HEADER_${headername}" var_headername)
  find_path(${var_headername} ${headername})
  if(NOT ${var_headername})
    message(FATAL_ERROR "Required Boost header ${headername} was not found.")
  else()
    include_directories(${${var_headername}})
  endif()
  mark_as_advanced(${var_headername})
endforeach()
