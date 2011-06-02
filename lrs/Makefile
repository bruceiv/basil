CC = gcc
CPPFLAGS = -DTIMES -DSIGNALS -DGMP -DLRS_QUIET
CXXFLAGS = -ggdb
LDFLAGS = -lgmpxx -lgmp


.PHONY:  clean

liblrs.a: lrs.cpp lrsMatrix.cpp lrs.hpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -o lrslib.o -c lrslib.c $(LDFLAGS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o lrsgmp.o -c lrsgmp.c $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o lrsMatrix.o -c lrsMatrix.cpp $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o lrs.o -c lrs.cpp $(LDFLAGS)
	$(AR) rcs liblrs.a lrs.o lrsMatrix.o lrslib.o lrsgmp.o

clean:  
	rm liblrs.a lrs.o lrsMatrix.o lrsgmp.o lrslib.o
