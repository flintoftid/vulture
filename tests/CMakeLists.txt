# 
# This file is part of Vulture.
#
# Vulture finite-difference time-domain electromagnetic solver.
# Copyright (C) 2011-2013 Ian David Flintoft
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
# Author: Ian Flintoft <ian.flintoft@york.ac.uk>
#

# Test driver function.
# Each test for this driver should have a subdirectory called TESTNAME
# containing a mesh called TESTNAME.mesh. Optiionally it can have:
#
# 1. Gnuplot script plottd.gnp.cmake to plot the time.
# 2. Gnuplot script plotfd.gnp.cmake to plot the frequency domain data.
# 3. ASCII field output files eh_<tag>_td.asc and eh_<tag>_fd.asc to
#    validate the outputs of the sover.
# 4. A cmake script process.cmake to create process.dat files and run thr
#    required processing tools.
#

function( vulture_test TESTNAME ) 

  configure_file( "${VULTURE_SOURCE_DIR}/tests/startup.m.cmake" startup.m )

  file( COPY ${TESTNAME}.mesh DESTINATION . )

  if( CHECK_LIMITS )
    add_test( NAME ${TESTNAME}  COMMAND ${VULTURE_BINARY_DIR}/src/vulture -v ${TESTNAME}.mesh )
  else( CHECK_LIMITS )
    add_test( NAME ${TESTNAME}_vulture  COMMAND ${VULTURE_BINARY_DIR}/src/vulture -v ${TESTNAME}.mesh )
    add_test( NAME ${TESTNAME}_gvulture COMMAND ${VULTURE_BINARY_DIR}/src/gvulture -p ${TESTNAME}.mesh )
    if( PROCESSING_TESTS )
      if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/process.cmake )
        configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/process.cmake
                        ${CMAKE_CURRENT_BINARY_DIR}/process.cmake )
        add_test( NAME ${TESTNAME}_process COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/process.cmake )
      endif()
    endif( PROCESSING_TESTS )
    if( PYTHON_TESTS )
      file( GLOB TD_DATAFILES "${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/eh_*_td.asc" "${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/wf_*_td.asc" )
      foreach( TD_DATAFILE ${TD_DATAFILES} )
        get_filename_component( basename ${TD_DATAFILE} NAME )
        string( LENGTH ${basename} namelen )
        math( EXPR lastchar "${namelen} - 10" )
        string( SUBSTRING ${basename} 3 ${lastchar} tag )
        add_test( NAME ${TESTNAME}_difftd_${tag} COMMAND ${PYTHON_EXECUTABLE} ${VULTURE_SOURCE_DIR}/util/pydiff.py ${TD_DATAFILE} ${basename} )
      endforeach()
      file( GLOB FD_DATAFILES "${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/eh_*_fd.asc" "${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/wf_*_fd.asc" )
      foreach( FD_DATAFILE ${FD_DATAFILES} )
        get_filename_component( basename ${FD_DATAFILE} NAME )
        string( LENGTH ${basename} namelen )
        math( EXPR lastchar "${namelen} - 10" )
        string( SUBSTRING ${basename} 3 ${lastchar} tag )
        #add_test( NAME ${TESTNAME}_difffd_${tag} COMMAND ${PYTHON_EXECUTABLE} ${VULTURE_SOURCE_DIR}/util/pydiff.py ${FD_DATAFILE} ${basename} )
      endforeach()
    endif( PYTHON_TESTS )
    if( OCTAVE_TESTS )
      if( EXISTS ${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/startup.m.cmake )
        configure_file( ${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/startup.m.cmake
                        ${VULTURE_BINARY_DIR}/tests/${TESTNAME}/startup.m )
      endif()
      file( GLOB MFILES "${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/*.m" )
      foreach( MFILE ${MFILES} )
        get_filename_component( tag ${MFILE} NAME_WE )
        add_test( NAME ${TESTNAME}_octave_${tag} COMMAND ${OCTAVE} ${MFILE} )
      endforeach()
    endif( OCTAVE_TESTS )
    if( GNUPLOT_TESTS )
      if( EXISTS ${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/plottd.gnp.cmake )
        configure_file( ${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/plottd.gnp.cmake
                        ${VULTURE_BINARY_DIR}/tests/${TESTNAME}/plottd.gnp )
        add_test( NAME ${TESTNAME}_gnuplot_td   COMMAND ${GNUPLOT_EXECUTABLE} plottd.gnp )
      endif()
      if( EXISTS ${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/plotfd.gnp.cmake )
        configure_file( ${VULTURE_SOURCE_DIR}/tests/${TESTNAME}/plotfd.gnp.cmake
                        ${VULTURE_BINARY_DIR}/tests/${TESTNAME}/plotfd.gnp )
        add_test( NAME ${TESTNAME}_gnuplot_fd   COMMAND ${GNUPLOT_EXECUTABLE} plotfd.gnp )
      endif()
      add_test( NAME ${TESTNAME}_gnuplot_mesh COMMAND ${GNUPLOT_EXECUTABLE} mesh.gnp )
    endif( GNUPLOT_TESTS )
  endif( CHECK_LIMITS )

endfunction()

# Test driver function for parser tests.
function( vulture_parser_test TESTNAME WILL_FAIL ) 

  file( COPY ${TESTNAME}.mesh DESTINATION . )

  add_test( NAME ${TESTNAME}  COMMAND ${VULTURE_BINARY_DIR}/src/vulture -v ${TESTNAME}.mesh )
  set_property( TEST ${TESTNAME} PROPERTY WILL_FAIL ${WILL_FAIL} )

endfunction()

#
# Parser tests.
#
add_subdirectory( parser )

#
# Basic propagation tests.
#
add_subdirectory( freespace_mur )
add_subdirectory( freespace_pml )

#
# Source tests.
#
add_subdirectory( softsrc_idx )
add_subdirectory( softsrc_idy )
add_subdirectory( softsrc_idz )
add_subdirectory( softsrc_imdx )
add_subdirectory( softsrc_imdy )
add_subdirectory( softsrc_imdz )
add_subdirectory( softsrc_sum )
add_subdirectory( hardsrc_ex )
add_subdirectory( hardsrc_ey )
add_subdirectory( hardsrc_ez )
add_subdirectory( hardsrc_sum )
add_subdirectory( hardsrc_jmsxy )

#
# Dipole tests.
#
add_subdirectory( dipole_vrx )
add_subdirectory( dipole_vry )
add_subdirectory( dipole_vrz )
add_subdirectory( dipole_sum )

#
# ABC tests.
#
add_subdirectory( sgref_mur )
add_subdirectory( sgref_pml6 )
add_subdirectory( sgref_pml12 )
if ( BIG_TESTS )
  add_subdirectory( dgref_ref )
  add_subdirectory( dgref_mur )
  add_subdirectory( dgref_pml6 )
  add_subdirectory( dgref_pml12 )
  add_subdirectory( dgrefmb_ref )
  add_subdirectory( dgrefmb_pml6 )
endif( BIG_TESTS )

#
# Test PMC/PEC boundaries.
#
add_subdirectory( parplate_empty_yzx )
add_subdirectory( parplate_empty_yxz )
add_subdirectory( parplate_empty_zxy )
add_subdirectory( parplate_empty_zyx )
add_subdirectory( parplate_empty_xyz )
add_subdirectory( parplate_empty_xzy )
add_subdirectory( parplate_empty_sum )

#
# Test PMC/PEC boundaries and 1D simulation.
# 
add_subdirectory( parplate1D_empty_yzx )
add_subdirectory( parplate1D_empty_yxz )
add_subdirectory( parplate1D_empty_zxy )
add_subdirectory( parplate1D_empty_zyx )
add_subdirectory( parplate1D_empty_xyz )
add_subdirectory( parplate1D_empty_xzy )
add_subdirectory( parplate1D_empty_sum )

#
# Basic simple medium test.
#
add_subdirectory( parplate1D_diel )
add_subdirectory( parplate1D_simpleslab_yzx )
add_subdirectory( parplate1D_simpleslab_yxz )
add_subdirectory( parplate1D_simpleslab_zxy )
add_subdirectory( parplate1D_simpleslab_zyx )
add_subdirectory( parplate1D_simpleslab_xyz )
add_subdirectory( parplate1D_simpleslab_xzy )
add_subdirectory( parplate1D_simpleslab_sum )

#
# PML termination of simple medium test.
#
add_subdirectory( parplate1D_simplepml_yzx )
add_subdirectory( parplate1D_simplepml_yxz )
add_subdirectory( parplate1D_simplepml_zxy )
add_subdirectory( parplate1D_simplepml_zyx )
add_subdirectory( parplate1D_simplepml_xyz )
add_subdirectory( parplate1D_simplepml_xzy )
add_subdirectory( parplate1D_simplepml_sum )

#
# Basic Debye medium test.
#
add_subdirectory( parplate1D_debyeslab_yzx )
add_subdirectory( parplate1D_debyeslab_yxz )
add_subdirectory( parplate1D_debyeslab_zxy )
add_subdirectory( parplate1D_debyeslab_zyx )
add_subdirectory( parplate1D_debyeslab_xyz )
add_subdirectory( parplate1D_debyeslab_xzy )
add_subdirectory( parplate1D_debyeslab_sum )

#
# PML termination of Debye medium.
#
add_subdirectory( parplate1D_debyepml_yzx )
add_subdirectory( parplate1D_debyepml_yxz )
add_subdirectory( parplate1D_debyepml_zxy )
add_subdirectory( parplate1D_debyepml_zyx )
add_subdirectory( parplate1D_debyepml_xyz )
add_subdirectory( parplate1D_debyepml_xzy )
add_subdirectory( parplate1D_debyepml_sum )
add_subdirectory( parplate1D_debyepml_ext )

#
# Basic planewave test.
#
add_subdirectory( planewave_empty_xyz )
add_subdirectory( planewave_empty_xzy )
add_subdirectory( planewave_empty_yxz )
add_subdirectory( planewave_empty_yzx )
add_subdirectory( planewave_empty_zyx )
add_subdirectory( planewave_empty_zxy )
add_subdirectory( planewave_empty_sum )

#
# Test partial planewave on PMC boundaries.
#
add_subdirectory( parplate1D_planewave_xyz )
add_subdirectory( parplate1D_planewave_xzy )
add_subdirectory( parplate1D_planewave_yxz )
add_subdirectory( parplate1D_planewave_yzx )
add_subdirectory( parplate1D_planewave_zyx )
add_subdirectory( parplate1D_planewave_zxy )
add_subdirectory( parplate1D_planewave_sum )

#
# Planewave cut by PEC.
#
add_subdirectory( planewave_pec_xyz )
add_subdirectory( planewave_pec_xzy )
add_subdirectory( planewave_pec_yxz )
add_subdirectory( planewave_pec_yzx )
add_subdirectory( planewave_pec_zyx )
add_subdirectory( planewave_pec_zxy )
add_subdirectory( planewave_pec_sum )

#
# Waveform tests.
#
add_subdirectory( waveform_ext )

#
# Internal surface tests.
#
add_subdirectory( closedbox_pec )

# SIBC tests
if( WITH_SIBC )

  add_subdirectory( parplate1D_sibcspar_xyz )
  add_subdirectory( parplate1D_sibcspar_xzy )
  add_subdirectory( parplate1D_sibcspar_yxz )
  add_subdirectory( parplate1D_sibcspar_yzx )
  add_subdirectory( parplate1D_sibcspar_zxy )
  add_subdirectory( parplate1D_sibcspar_zyx )
  add_subdirectory( parplate1D_sibcspar_sum )
  add_subdirectory( parplate1D_sibcspar_ext )

  add_subdirectory( parplate1D_sibcgo04_xyz )
  add_subdirectory( parplate1D_sibcgo04_xzy )
  add_subdirectory( parplate1D_sibcgo04_yxz )
  add_subdirectory( parplate1D_sibcgo04_yzx )
  add_subdirectory( parplate1D_sibcgo04_zxy )
  add_subdirectory( parplate1D_sibcgo04_zyx )
  add_subdirectory( parplate1D_sibcgo04_sum )

  add_subdirectory( closedbox_spar )
  add_subdirectory( closedbox_go04 )

endif( WITH_SIBC )

