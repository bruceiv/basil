#   Copyright: Aaron Moss, 2012, moss.aaron@unb.ca

#   This file is part of Basil.
#
#   Basil is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as
#   published by the Free Software Foundation, either version 3 of
#   the License, or (at your option) any later version.
#
#   Basil is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with Basil.  If not, see <http://www.gnu.org/licenses/>.

CPPFLAGS = -DTIMES -DGMP -DLRS_QUIET -DLRS_THREADSAFE -DBAS_WALLTIME
CXXFLAGS = -O2 -Wall -Wno-unused -fopenmp -I .
#CXXFLAGS = -ggdb -O0 -Wall -Wno-unused -fopenmp -I .
LDFLAGS = -lboost_program_options -Llrs -llrs -lgmpxx -lgmp

# object files to include in this executable
OBJS = automorphism.o parse.o gram.o dfs.o metric.o fund_domain.o
# object files for parallel version of executable
OBJSP = automorphism.o parse.o gram.o dfsp.o metric.o fund_domain.o

# rules for constructions of objects from sources
.cpp.o:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(LDFLAGS)

.PHONY:  clean clean_p clean_all clean_doc doc lrs

# generate main program
basil:  lrs $(OBJS) main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o basil main.cpp $(OBJS) $(LDFLAGS)

# generate multithreaded main program
dfsp.o:  dfsp.cpp
	$(CXX) $(CPPFLAGS) -DBAS_MT $(CXXFLAGS)	-c dfsp.cpp $(LDFLAGS)

basilp:  lrs $(OBJSP) main.cpp
	$(CXX) $(CPPFLAGS) -DBAS_MT $(CXXFLAGS)	-o basilp main.cpp $(OBJSP) $(LDFLAGS)

# generate lrs library
lrs:  
	cd lrs && make

# clean generated files
clean:  
	-rm $(OBJS) basil

clean_p:
	-rm $(OBJSP) basilp

# clean all generated files (including libraries and documentation)
clean_all:  clean clean_p clean_doc
	-cd lrs && make clean

#clean documentation
clean_doc:  
	rm -rf doc/*

# generate documentation (requires doxygen)
doc:  
	doxygen Doxyfile

# package code
tarball:
	cd .. && tar -zcvf basil.tar.gz --exclude basil/.git basil/
