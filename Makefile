LDFLAGS = -Llrs -lgmpxx -lgmp -llrs


.PHONY:  clean clean_all clean_doc doc lrs

# generate main program
basil:  main.cpp lrs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o basil main.cpp $(LDFLAGS)

# generate lrs library
lrs:  
	cd lrs && make

# clean generated files
clean:  
	rm basil

# clean all generated files (including libraries and documentation)
clean_all:  clean clean_doc
	cd lrs && make clean

#clean documentation
clean_doc:  
	rm -rf doc/*

# generate documentation (requires doxygen)
doc:  
	doxygen Doxyfile