CPPFLAGS = -DTIMES -DGMP -DLRS_QUIET -DLRS_THREADSAFE
CXXFLAGS = -ggdb -O2 -Wall -Wno-unused
LDFLAGS = -lboost_program_options-mt -Llrs -llrs -lgmpxx -lgmp

# object files to include in this executable
OBJS = automorphism.o parse.o gram.o dfs.o
# object files for parallel version of executable
OBJSP = automorphism.o parse.o gram.o dfsp.o

# rules for constructions of objects from sources
.cpp.o:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(LDFLAGS)

.PHONY:  clean clean_all clean_doc doc lrs

# generate main program
basil:  lrs $(OBJS) main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o basil main.cpp $(OBJS) $(LDFLAGS)

# generate multithreaded main program
dfsp.o:  dfsp.cpp
	$(CXX) $(CPPFLAGS) -DBAS_MT $(CXXFLAGS) -fopenmp \
	-c dfsp.cpp $(LDFLAGS)

basilp:  lrs $(OBJSP) main.cpp
	$(CXX) $(CPPFLAGS) -DBAS_MT $(CXXFLAGS) -fopenmp \
	-o basilp main.cpp $(OBJSP) $(LDFLAGS)

# generate lrs library
lrs:  
	cd lrs && make

# clean generated files
clean:  
	-rm $(OBJS) basil basilp

# clean all generated files (including libraries and documentation)
clean_all:  clean clean_doc
	-cd lrs && make clean

#clean documentation
clean_doc:  
	rm -rf doc/*

# generate documentation (requires doxygen)
doc:  
	doxygen Doxyfile
