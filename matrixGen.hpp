#ifndef _MATRIX_GEN_HPP_
#define _MATRIX_GEN_HPP_

#include <istream>

#include <gmpxx.h>

#include "basilCommon.hpp"

namespace basil {
	
	/** Allocates new matrix on heap and returns it.
	 *  expects whitespace-delimited values as follows: first m and n, the 
	 *  dimensions of the matrix, then the m*n data values of the matrix.
	 */
	matrix_ptr genMatrixFromStream(std::istream& in) {
		ind n, d;
		in >> n;
		in >> d;
		
		matrix_ptr m(new matrix(n, d));

		mpq_class t;
		for (ind i = 0; i < n; i++) {
			for (ind j = 0; j < d; j++) {
				in >> t;
				(*m)[i][j] = t;
			}
		}
		
		return m;
	}
	
} /* namespace basil */

#endif /* _MATRIX_GEN_HPP_ */