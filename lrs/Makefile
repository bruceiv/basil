# Copyright: Aaron Moss, 2012, moss.aaron@unb.ca

# This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


CC = gcc
CPPFLAGS = -DTIMES -DGMP -DLRS_QUIET -DLRS_THREADSAFE
CFLAGS = -O2 -fopenmp
#CFLAGS = -ggdb -O0 -fopenmp
CXXFLAGS = $(CFLAGS) -Wall -Wno-unused
LDFLAGS = -lgmpxx -lgmp

# object files to include in this library
OBJS = lrslib.o lrsgmp.o matrix.o lrs.o

# rules for constructions of objects from sources
.c.o:  
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< $(LDFLAGS)
.cpp.o:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(LDFLAGS)

.PHONY:  clean

liblrs.a: $(OBJS)
	$(AR) rcs liblrs.a $(OBJS)

clean:  
	-rm $(OBJS) liblrs.a
