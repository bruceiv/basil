/** Implements lrs::matrix class from matrix.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <istream>
#include <new>
#include <ostream>

#include <gmpxx.h>

#include "matrix.hpp"


namespace lrs {
	
	////////////////////////////////////////////////////////////////////////////
	// Constructors, destructor, and assignment operator
	////////////////////////////////////////////////////////////////////////////
	
	matrix::matrix(ind n, ind d) throw(std::bad_alloc)
			: n_(n), d_(d), 
			num(new matrix_t(lrs_alloc_mp_matrix(n, d))), 
			den(new matrix_t(lrs_alloc_mp_matrix(n, d))) { 
		if (*num == 0 || *den == 0) throw std::bad_alloc();
	}
	
	matrix::matrix(matrix const& that) throw(std::bad_alloc)
			: n_(that.n_), d_(that.d_), 
			num(new matrix_t(lrs_alloc_mp_matrix(that.n_, that.d_))), 
			den(new matrix_t(lrs_alloc_mp_matrix(that.n_, that.d_))) { 
		
		if (*num == 0 || *den == 0) throw std::bad_alloc();
		
		for (ind i = 0; i < n_; i++) {
			for (ind j = 0; j < d_; j++) {
				mpz_set((*num)[i][j], (*that.num)[i][j]);
				mpz_set((*den)[i][j], (*that.den)[i][j]);
			}
		}
	}
	
	matrix::~matrix() {
		lrs_clear_mp_matrix(*num, n_, d_);
		lrs_clear_mp_matrix(*den, n_, d_);
		delete num;
		delete den;
	}
	
	matrix& matrix::operator= (matrix const& that) throw(std::bad_alloc) {
		//rebuild matrix storage if assigning a matrix of another dimension
		if (n_ != that.n_ || d_ != that.d_) {
			lrs_clear_mp_matrix(*num, n_, d_);
			lrs_clear_mp_matrix(*den, n_, d_);
			delete num;
			delete den;
			
			n_ = that.n_; d_ = that.d_;
			num = new matrix_t(lrs_alloc_mp_matrix(n_, d_));
			den = new matrix_t(lrs_alloc_mp_matrix(n_, d_));
			
			if (*num == 0 || *den == 0) throw std::bad_alloc();
		}
		
		for (ind i = 0; i < n_; i++) {
			for (ind j = 0; j < d_; j++) {
				mpz_set((*num)[i][j], (*that.num)[i][j]);
				mpz_set((*den)[i][j], (*that.den)[i][j]);
			}
		}
		
		return *this;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Matrix operators
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< (std::ostream& out, matrix const & m) {
		//print dimensions
		out << "(" << m.n_ << "," << m.d_ << ")";
		
		//print data values
		bool firstRow = true;
		out << "[";
		for (ind i = 0; i < m.n_; i++) {
			if (!firstRow) out << ","; else firstRow = false;
			bool firstCol = true;
			out << "[";
			for (ind j = 0; j < m.d_; j++) {
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
	
	/* inline */ vector_t& matrix::matrix_row::num() {
		return (*m->num)[i];
	}
	
	/* inline */ vector_t& matrix::matrix_row::den() {
		return (*m->den)[i];
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Indirection operators
	////////////////////////////////////////////////////////////////////////////
	
	/* inline */ matrix::matrix_ind::operator mpq_class () const {
		return mpq_class(mpz_class((*m->num)[i][j]), 
						 mpz_class((*m->den)[i][j]));
	}
	
	/* inline */ matrix::const_matrix_ind::operator mpq_class () const {
		return mpq_class(mpz_class((*m->num)[i][j]), 
						 mpz_class((*m->den)[i][j]));
	}
	
	/* inline */ matrix::matrix_ind& matrix::matrix_ind::operator= (
			const mpq_class& that) {
		mpz_set((*m->num)[i][j], that.get_num_mpz_t());
		mpz_set((*m->den)[i][j], that.get_den_mpz_t());
		return *this;
	}
	
	/* inline */ std::istream& operator>> (std::istream& in, 
									 matrix::matrix_ind& m_i) {
		mpq_class t; in >> t; m_i = t;
	}
	
	/* inline */ std::ostream& operator<< (std::ostream& out,
									 matrix::matrix_ind& m_i) {
		out << mpq_class(m_i);
	}
	
	/* inline */ std::ostream& operator<< (std::ostream& out,
									 matrix::const_matrix_ind const & m_i) {
		out << mpq_class(m_i);
	}
	
} /* namespace lrs */
