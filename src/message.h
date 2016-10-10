/* 
 * This file is part of Vulture.
 *
 * Vulture finite-difference time-domain electromagnetic solver.
 * Copyright (C) 2011-2016 Ian David Flintoft
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 * Author: Ian Flintoft <ian.flintoft@googlemail.com>
 *
 */

#ifndef _MESSAGE_H_
#define _MESSAGE_H_

#include <errno.h>

typedef enum msg_type {

  MSG_DEBUG3,
  MSG_DEBUG2,
  MSG_DEBUG1,
  MSG_INFO,
  MSG_LOG,
  MSG_WARN,
  MSG_ERROR

} MessageType;

void startMessaging( char *logFileName , MessageType minimumLevel , char *progName , int versionMajor , int versionMinor , int versionPatch );
void message(  MessageType status , int errnum , const char *message , ... );
void stopMessaging( void );

#endif


