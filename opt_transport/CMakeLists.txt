# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( opt_transport )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package


include_directories ("${PROJECT_SOURCE_DIR}/measures")
add_subdirectory (measures)


# Creating entries for target: opt_transport
# ############################

add_executable( opt_transport  opt_transport.cpp )

# target_compile_features(opt_transport PRIVATE cxx_range_for)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS opt_transport )

# Link the executable to CGAL and third-party libraries
target_link_libraries(opt_transport   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} lbfgs Source Target )

install (TARGETS opt_transport DESTINATION ${PROJECT_SOURCE_DIR}/bin)
