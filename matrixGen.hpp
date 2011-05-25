#ifndef _MATRIX_GEN_HPP_
#define _MATRIX_GEN_HPP_

#include <istream>

#include "basilCommon.hpp"

namespace basil {
	
	//allocates new matrix on heap and returns it.
	//expects whitespace-delimited values as follows: first n and m, the 
	// dimensions of the matrix (ind), then the n*m data values of the matrix.
	shared_ptr<matrix> genMatrixfromStream(std::istream& in) {
		ind n, d;
		in >> n;
		in >> d;
		
		shared_ptr<matrix> m(new matrix(n, d));
			
		for (ind i = 0; i < n; i++) {
			for (ind j = 0; j < d; j++) {
				in >> (*m) (i, j);
			}
		}
		
		return m;
	}
	
}

#endif /* _MATRIX_GEN_HPP_ */