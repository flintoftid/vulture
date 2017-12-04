# 
# This file is part of Vulture.
#
# Vulture finite-difference time-domain electromagnetic solver.
# Copyright (C) 2011-2016 Ian David Flintoft
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#
# Author: Ian Flintoft <ian.flintoft@googlemail.com>
#

project ( VULTURE C )

# OpenMP support flaky before 2.8.6?
cmake_minimum_required( VERSION 2.8.6 )

# Add path to local CMake modules.
set( CMAKE_MODULE_PATH "${VULTURE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

# Project name and version.
set( ${PROJECT_NAME}_PACKAGE        "VULTURE" )
set( ${PROJECT_NAME}_VERSION_MAJOR  "0" )
set( ${PROJECT_NAME}_VERSION_MINOR  "7" )
set( ${PROJECT_NAME}_VERSION_PATCH  "1" )
set( ${PROJECT_NAME}_VERSION_STRING "\"${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}\"" )

# Build options.
option( CHECK_LIMITS       "Build limit checking version"   OFF )
option( WITH_OPENMP        "Compile with OpenMP support"    OFF )
option( WITH_SIBC          "Compile with SIBC support"      OFF )
option( USE_INDEXED_MEDIA  "Compile using indexed media"    OFF )
option( USE_SCALED_FIELDS  "Compile using scaled fields"    OFF )
option( USE_AVERAGED_MEDIA "Compile using media averaging"  OFF )

option( WITH_TESTS         "Configure test suite"           OFF )
option( GNUPLOT_TESTS      "Add gnuplot tests"              OFF )
option( PYTHON_TESTS       "Add python validation tests"    OFF )
option( OCTAVE_TESTS       "Add octave validation tests"    OFF )
option( PROCESSING_TESTS   "Add processing tool tests"      OFF )
option( BIG_TESTS          "Add big RAM/CPU time tests"     OFF )

option( WITH_LATEX         "Compile LaTeX user manual"      OFF )

#
# Compile options.
#

# No easy way to enable/find C99 compiler in CMake. 
if( CMAKE_COMPILER_IS_GNUCC )
  # gcc/mingw32.
  message( STATUS "Using GNU C compiler ${C_COMPILER_NAME} ..." )
  set( CMAKE_C_FLAGS         "-Wall -pedantic-errors -std=c99" )
  set( CMAKE_C_FLAGS_PROFILE "-O0 -g -fprofile-arcs -ftest-coverage" )
  set( CMAKE_C_FLAGS_DEBUG   "-DDEBUG -O1 -g -fbounds-check" )
  #set( CMAKE_C_FLAGS_RELEASE "-m64 -Ofast -flto -DNDEBUG -mfpmath=sse -funroll-loops -march=core2" )
  set( CMAKE_C_FLAGS_RELEASE " -Ofast -flto -DNDEBUG -mfpmath=sse -funroll-loops -march=core2" )
  set( CMAKE_AR "gcc-ar" )
  set( CMAKE_RANLIB "gcc-ranlib" )
  # Note if we are using MingW so we can workaround some bugs.
  if( MINGW )
    add_definitions( -DUSE_MINGW_FIX )
    message( STATUS "Applying MingW vprintf workaround" ) 
  endif( MINGW )
else( CMAKE_COMPILER_IS_GNUCC )
  # Generic
  message( STATUS "Using generic C ${C_COMPILER_NAME} ..." )
  set( CMAKE_C_FLAGS         "" )
  set( CMAKE_C_FLAGS_PROFILE "-O0" )
  set( CMAKE_C_FLAGS_DEBUG   "-DDEBUG -O1 -g" )
  set( CMAKE_C_FLAGS_RELEASE "-O3 -DNDDEBUG" )
endif( CMAKE_COMPILER_IS_GNUCC )

#
# Check build option consistency and dependencies.
#

if( CHECK_LIMITS )
  if( USE_AVERAGED_MEDIA )
    message( FATAL_ERROR "Averaged media not compatiable with check limits" ) 
  else( USE_AVERAGED_MEDIA )
    add_definitions( -DCHECK_LIMITS )
  endif( USE_AVERAGED_MEDIA)
endif( CHECK_LIMITS )

if( USE_INDEXED_MEDIA )
  if( USE_SCALED_FIELDS OR USE_AVERAGED_MEDIA )
    message( FATAL_ERROR "Indexed media not compatiable with scaled fields and averaged media" ) 
  else( USE_SCALED_FIELDS OR USE_AVERAGED_MEDIA )
    add_definitions( -DUSE_INDEXED_MEDIA )
  endif( USE_SCALED_FIELDS OR USE_AVERAGED_MEDIA )
endif( USE_INDEXED_MEDIA )
 
if( USE_SCALED_FIELDS )
  if( USE_INDEXED_MEDIA )
    message( FATAL_ERROR "Scaled fields not compatiable with indexed media" ) 
  else( USE_INDEXED_MEDIA )
    add_definitions( -DUSE_SCALED_FIELDS )
  endif( USE_INDEXED_MEDIA )
endif( USE_SCALED_FIELDS )

if( USE_AVERAGED_MEDIA )
  if( USE_INDEXED_MEDIA )
    message( FATAL_ERROR "Averaged media not compatiable with indexed media" ) 
  else( USE_INDEXED_MEDIA )
    add_definitions( -DUSE_AVERAGED_MEDIA )
  endif( USE_INDEXED_MEDIA )
endif( USE_AVERAGED_MEDIA )

if( WITH_OPENMP )
  find_package( OpenMP ) 
  if( OPENMP_FOUND )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
    add_definitions( -DWITH_OPENMP )
  endif( OPENMP_FOUND )
endif( WITH_OPENMP )

if( WITH_SIBC )
  add_definitions( -DWITH_SIBC )
endif( WITH_SIBC )

# 
# Main souce code.
#

include_directories( uthash )
add_subdirectory ( src )

install( DIRECTORY m/ DESTINATION share/octave/packages/vulture FILES_MATCHING PATTERN "*.m" )

#
# Documentation.
#

install( FILES "${VULTURE_SOURCE_DIR}/Licence.txt"
               "${VULTURE_SOURCE_DIR}/ReadMe.md"
         DESTINATION share/doc/vulture
         COMPONENT instructionData )

install( DIRECTORY "${VULTURE_SOURCE_DIR}/examples/"
         DESTINATION share/doc/vulture/examples )

add_subdirectory ( doc )

#
# Testing.
#

if( WITH_TESTS )

  include( CTest )
  enable_testing()

  # Check if for Python.
  if( PYTHON_TESTS )
    find_package( PythonInterp )
    if( PYTHONINTERP_FOUND )
      message( STATUS "Found python interpreter: ${PYTHON_EXECUTABLE}" )
    else( PYTHONINTERP_FOUND )
      message( WARNING "Python interpreter not found: disabling diffs in tests" )
      option( PYTHON_TESTS OFF )
    endif( PYTHONINTERP_FOUND )
  endif( PYTHON_TESTS )

  # If running matlab tests requested:
  if( OCTAVE_TESTS )
    find_program( OCTAVE octave )
    if( OCTAVE )
      message( STATUS "Octave: ${OCTAVE}" )
      set( MEXEC ${OCTAVE} "-q" "--eval" CACHE INTERNAL "" )
    else( OCTAVE )
      find_program( MATLAB matlab )
      if( MATLAB )
        message( STATUS "MATLAB: ${MATALB}" )
        set( MEXEC  matlab "-nosplash" "-nojvm" "-nodisplay" "-nodesktop" "-r" CACHE INTERNAL "" ) 
      else( MATLAB )
        message( FATAL_ERROR "Cannot find Octave or MATLAB" )
      endif( MATLAB )
    endif( OCTAVE )
  endif( OCTAVE_TESTS )

  # Check for gnuplot.
  if( GNUPLOT_TESTS )
    find_package( Gnuplot )
    if( GNUPLOT_FOUND )
      message( STATUS "Found gnuplot: ${GNUPLOT_EXECUTABLE}" )
    else( GNUPLOT_FOUND )
      message( WARNING "Gnuplot not found: disabling gnuplot in tests" )
      option( GNUPLOT_TESTS OFF )
    endif( GNUPLOT_FOUND )
  endif( GNUPLOT_TESTS )

  # Check for xtime, xtransall, xfreq
  if( PROCESSING_TESTS )
    find_program( XTIME_EXECUTABLE NAMES xtime )
    find_program( XTRANSALL_EXECUTABLE NAMES xtransall )
    find_program( XFREQ_EXECUTABLE NAMES xfreq )
    if( XTIME_EXECUTABLE AND XTRANSALL_EXECUTABLE AND  XFREQ_EXECUTABLE )
      message( STATUS "Found processing tools (xtime): ${XTIME_EXECUTABLE}" )  
      message( STATUS "Found processing tools (xtransall): ${XTRANSALL_EXECUTABLE}" )  
      message( STATUS "Found processing tools (xfreq): ${XFREQ_EXECUTABLE}" )  
    else()
      message( WARNING "Processing tools not all found: disabling processing tests" )
      option( PROCESSING_TESTS OFF )  
    endif()
  endif( PROCESSING_TESTS )

  add_subdirectory ( tests )

endif( WITH_TESTS )

#
# Package generation.
#

set( CPACK_PACKAGE_NAME                "vulture" )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Vulture FDTD code" )
set( CPACK_PACKAGE_VENDOR              "Ian Flintoft <ian.flintoft@googlemail.com>" )
set( CPACK_PACKAGE_CONTACT             "Ian Flintoft <ian.flintoft@googlemail.com>" )
set( CPACK_PACKAGE_DESCRIPTION_FILE    "${VULTURE_SOURCE_DIR}/ReadMe.md" )
set( CPACK_RESOURCE_FILE_LICENSE       "${VULTURE_SOURCE_DIR}/Licence.txt" )
set( CPACK_PACKAGE_VERSION_MAJOR       "${VULTURE_VERSION_MAJOR}" )
set( CPACK_PACKAGE_VERSION_MINOR       "${VULTURE_VERSION_MINOR}" )
set( CPACK_PACKAGE_VERSION_PATCH       "${VULTURE_VERSION_PATCH}" )
set( CPACK_PACKAGE_VERSION             "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )
set( CPACK_PACKAGE_INSTALL_DIRECTORY   "${CPACK_PACKAGE_NAME}-${CPACK_VERSION}" )
set( CPACK_SOURCE_IGNORE_FILES         "/\\\\.hg/;\\\\.swp$;\\\\.#;/#;\\\\.*~;\\\\.hgtags" )
set( CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR}" )

if( UNIX )  
   set( CPACK_GENERATOR "TGZ" )
   set( CPACK_SOURCE_GENERATOR "TGZ" )
elseif( WIN32 AND NOT UNIX )
   set( CPACK_GENERATOR "NSIS" )
   set( CPACK_SOURCE_GENERATOR "ZIP" )
endif( UNIX)

include( CPack )
