CPPFLAGS = -DTIMES -DGMP -DLRS_QUIET
CXXFLAGS = -ggdb -O2 -Wall -Wno-unused
LDFLAGS = -lboost_program_options -Llrs -llrs -lgmpxx -lgmp

# object files to include in this executable
OBJS = parse.o gram.o dfs.o

# rules for constructions of objects from sources
.cpp.o:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(LDFLAGS)


.PHONY:  clean clean_all clean_doc doc lrs

# generate main program
basil:  lrs $(OBJS) main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o basil main.cpp $(OBJS) $(LDFLAGS)

# generate lrs library
lrs:  
	cd lrs && make

# clean generated files
clean:  
	-rm $(OBJS) basil

# clean all generated files (including libraries and documentation)
clean_all:  clean clean_doc
	-cd lrs && make clean

#clean documentation
clean_doc:  
	rm -rf doc/*

# generate documentation (requires doxygen)
doc:  
	doxygen Doxyfile