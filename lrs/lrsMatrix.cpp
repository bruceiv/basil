/** Implements lrs::matrix class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <istream>
#include <ostream>

#include <gmpxx.h>

#include "lrs.hpp"
extern "C" {
#include LRSXX_MP_H
}

namespace lrs {
	
	////////////////////////////////////////////////////////////////////////////
	// Constructors and destructors
	////////////////////////////////////////////////////////////////////////////
	
	matrix::matrix(ind n, ind d) 
		: n(n), d(d), 
		num(lrs_alloc_mp_matrix(n, d)), den(lrs_alloc_mp_matrix(n, d)) { }
	
	matrix::~matrix() {
		lrs_clear_mp_matrix(num, n, d);
		lrs_clear_mp_matrix(den, n, d);
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Matrix operators
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< (std::ostream& out, matrix& m) {
		//print dimensions
		out << "(" << m.n << "," << m.d << ")";
		
		//print data values
		bool firstRow = true;
		out << "[";
		for (ind i = 0; i < m.n; i++) {
			if (!firstRow) out << ","; else firstRow = false;
			bool firstCol = true;
			out << "[";
			for (ind j = 0; j < m.d; j++) {
				if (!firstCol) out << ","; else firstCol = false;
				out << m[i][j];
			}
			out << "]";
		}
		out << "]";
		
		return out;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Matrix row operators
	////////////////////////////////////////////////////////////////////////////
	
	inline vector_t& matrix::matrix_row::num() {
		return m->num[i];
	}
	
	inline vector_t& matrix::matrix_row::den() {
		return m->den[i];
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Indirection operators
	////////////////////////////////////////////////////////////////////////////
	
	inline matrix::matrix_ind::operator mpq_class () const {
		return mpq_class(mpz_class(m->num[i][j]), mpz_class(m->den[i][j]));
	}
	
	inline matrix::const_matrix_ind::operator mpq_class () const {
		return mpq_class(mpz_class(m->num[i][j]), mpz_class(m->den[i][j]));
	}
	
	inline matrix::matrix_ind& matrix::matrix_ind::operator= (
			const mpq_class& that) {
		mpz_set(m->num[i][j], that.get_num_mpz_t());
		mpz_set(m->den[i][j], that.get_den_mpz_t());
		return *this;
	}
	
	inline std::istream& operator>> (std::istream& in, 
									 matrix::matrix_ind& m_i) {
		mpq_class t; in >> t; m_i = t;
	}
	
	inline std::ostream& operator<< (std::ostream& out,
									 matrix::matrix_ind& m_i) {
		out << mpq_class(m_i);
	}
	
	inline std::ostream& operator<< (std::ostream& out,
									 matrix::const_matrix_ind& m_i) {
		out << mpq_class(m_i);
	}
	
} /* namespace lrs */
