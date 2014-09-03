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

if( WITH_SIBC )
  set( SIBC_SOURCES sibc.c filter.c )
  set( SIBC_INCLUDES sibc.h filter.h )
else( WITH_SIBC)
  set( SIBC_SOURCES )
  set( SIBC_INCLUDES )
endif( WITH_SIBC )

set( VULTURE_SOURCES  fdtd_types.c physical.c message.c alloc_array.c simulation.c   
                      bounding_box.c mesh.c grid.c pml.c gnuplot.c gmsh.c timer.c memory.c
                      medium.c block.c boundary.c surface.c waveform.c source.c planewave.c 
                      observer.c util.c mur.c debye.c wire.c line.c ${SIBC_SOURCES} )

set( VULTURE_INCLUDES fdtd_types.h physical.h message.h alloc_array.h simulation.h
                      bounding_box.h mesh.h grid.h pml.h gnuplot.h gmsh.h vulture.h timer.h memory.h
                      medium.h block.h boundary.h surface.h waveform.h source.h planewave.h 
                      observer.h util.h mur.h debye.h wire.h line.h ${SIBC_INCLUDES} )

add_library( vult STATIC ${VULTURE_SOURCES} )
  
add_executable( vulture vulture.c )
target_link_libraries( vulture vult m )

add_executable( gvulture gvulture.c )
target_link_libraries( gvulture vult m )

install( TARGETS vulture gvulture RUNTIME DESTINATION bin )
