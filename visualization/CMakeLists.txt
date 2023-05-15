# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( create_diagram )


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


# Creating entries for target: create_diagram
# ############################

add_executable( create_diagram  create_diagram.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS create_diagram )

# Link the executable to CGAL and third-party libraries
target_link_libraries(create_diagram PRIVATE CGAL::CGAL )