CPPFLAGS = 
CXXFLAGS = -ggdb -Wall -Wno-unused
LDFLAGS = -lboost_unit_test_framework

# rules for constructions of objects from sources
.cpp.out:  
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

.PHONY:  clean test

test:  test.out
	@echo "*** LRU cache test ***"
	@./test.out

clean:  
	rm test.out
