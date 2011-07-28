CC = gcc
CPPFLAGS = -DTIMES -DSIGNALS -DGMP -DLRS_QUIET
CXXFLAGS = -O0 -ggdb -Wall -Wno-unused
LDFLAGS = -L../lrs -llrs -lboost_unit_test_framework -lgmpxx -lgmp

# rules for constructions of objects from sources
.cpp.out:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

.PHONY:  test clean


test:  clean test.out ../gram.hpp ../gram.cpp
	@echo "*** Gram matrix test ***"
	@./test.out

clean:  
	rm test.out