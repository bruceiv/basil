CC = gcc
CPPFLAGS = -DTIMES -DSIGNALS -DGMP -DLRS_QUIET
CXXFLAGS = -O0 -ggdb -Wall -Wno-unused
LDFLAGS = -L.. -llrs -lboost_unit_test_framework -lgmpxx -lgmp

# rules for constructions of objects from sources
.cpp.out:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

.PHONY:  test


test:  test.out ../matrix.hpp ../matrix.cpp
	@echo "*** LRS wrapper test ***"
	@./test.out
	@rm test.out
