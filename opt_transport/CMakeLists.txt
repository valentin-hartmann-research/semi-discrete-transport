# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( opt_transport )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


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

add_to_cached_list( CGAL_EXECUTABLE_TARGETS opt_transport )

# Link the executable to CGAL and third-party libraries
target_link_libraries(opt_transport PRIVATE CGAL::CGAL lbfgs Source Target)