CPPFLAGS = -DTIMES -DSIGNALS -DGMP -DLRS_QUIET
CXXFLAGS = -ggdb
LDFLAGS = -Llrs -llrs -lgmpxx -lgmp

# object files to include in this executable
OBJS = dfs.o

# rules for constructions of objects from sources
.cpp.o:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(LDFLAGS)


.PHONY:  clean clean_all clean_doc doc lrs

# generate main program
basil:  main.cpp $(OBJS) lrs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o basil main.cpp $(OBJS) $(LDFLAGS)

# generate lrs library
lrs:  
	cd lrs && make

# clean generated files
clean:  
	rm basil $(OBJS)

# clean all generated files (including libraries and documentation)
clean_all:  clean clean_doc
	cd lrs && make clean

#clean documentation
clean_doc:  
	rm -rf doc/*

# generate documentation (requires doxygen)
doc:  
	doxygen Doxyfile