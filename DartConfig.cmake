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

# Set site to hostname.
find_program( HOSTNAME_CMD NAMES hostname )
exec_program( ${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME )
set( SITE "${HOSTNAME}" )

# Use uname to construct build name.
find_program( UNAME NAMES uname )

macro( getuname name flag )
  exec_program( "${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}" )
endmacro( getuname )

getuname( osname -s )
getuname( osrel -r )
getuname( cpu -m )
set( BUILD_NAME "${osname}-${cpu}" )
set( CTEST_BUILD_NAME "${osname}-${cpu}" )

# Coverage and dynamic memory checking options.
find_program( CTEST_COVERAGE_COMMAND NAMES gcov )
find_program( CTEST_MEMORYCHECK_COMMAND NAMES valgrind )
#set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --tool=memcheck --leak-check=yes --show-reachable=yes --verbose --demangle=yes")
set( MEMORYCHECK_COMMAND_OPTIONS "--quiet --tool=memcheck --leak-check=yes --demangle=yes" )
