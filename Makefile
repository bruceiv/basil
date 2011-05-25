LDFLAGS = -Llrs -lgmpxx -lgmp -llrs


.PHONY:  clean clean_all lrs

basil:  main.cpp lrs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o basil main.cpp $(LDFLAGS)

lrs:  
	cd lrs && make

clean:  
	rm basil

clean_all:  clean
	cd lrs && make clean
