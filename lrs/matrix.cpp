/** Implements matrix and vector types from matrix.hpp.
 *
 *  @author Aaron Moss
 */

#include <iostream>
#include <cstdlib>
#include <istream>
#include <new>
#include <ostream>
#include <stdexcept>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include "cobasis.hpp"
#include "matrix.hpp"


namespace lrs {
	
	////////////////////////////////////////////////////////////////////////////
	//
	// VECTOR_MPQ_BASE CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpq_base::vector_mpq_base( mpq_class* v, ind d ) : v(v), d(d) {}
	
	////////////////////////////////////////////////////////////////////////////
	// Element access operators
	////////////////////////////////////////////////////////////////////////////
	
	mpq_class* vector_mpq_base::begin() { return v; }
	
	mpq_class const* vector_mpq_base::begin() const { return v; }
	
	mpq_class* vector_mpq_base::end() { return v + d; }
	
	mpq_class const* vector_mpq_base::end() const { return v + d; }
	
	ind vector_mpq_base::size() const { return d; }
	
	mpq_class& vector_mpq_base::operator[] ( ind i ) { return v[i]; }
	
	const mpq_class& vector_mpq_base::operator[] ( ind i ) const 
		{ return v[i]; }
	
	////////////////////////////////////////////////////////////////////////////
	// Output operator
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< ( std::ostream& o, vector_mpq_base const& v ) {
		o << "[";
		for (ind i = 0; i < v.d; i++) {
			o << " " << v.v[i];
		}
		o << " ]";
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Comparison operators
	////////////////////////////////////////////////////////////////////////////
	
	bool operator<  ( vector_mpq_base const& a, vector_mpq_base const& b ) 
		{ return compare(a,b) < 0; }
	bool operator== ( vector_mpq_base const& a, vector_mpq_base const& b ) 
		{ return compare(a,b) == 0; }
	bool operator>  ( vector_mpq_base const& a, vector_mpq_base const& b ) 
		{ return compare(a,b) > 0; }
	bool operator<= ( vector_mpq_base const& a, vector_mpq_base const& b ) 
		{ return compare(a,b) <= 0; }
	bool operator!= ( vector_mpq_base const& a, vector_mpq_base const& b ) 
		{ return compare(a,b) != 0; }
	bool operator>= ( vector_mpq_base const& a, vector_mpq_base const& b ) 
		{ return compare(a,b) >= 0; }
	
	int compare( vector_mpq_base const& a, vector_mpq_base const& b ) {
		mpq_class t; int s = 0;
		
		for (ind i = 0; i < a.d && i < b.d; ++i) {
			t = a.v[i] - b.v[i];
			s = sgn(t);
			if (s != 0) return s;
		}
		// if it reaches here, the two are lexicographically equal up to the 
		// end of the shorter vector (and s == 0)
		
		if (a.d < b.d) s = -1; else if (a.d > b.d) s = 1; /* else s = 0; */
		
		return s;
	}
	
	bool is_zero( vector_mpq_base const& v ) {
		for (ind i = 0; i < v.d; ++i) {
			if ( v.v[i] != 0 ) return false;
		}
		return true;
	}

	////////////////////////////////////////////////////////////////////////////
	// Mathematical operators
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpq operator+ (vector_mpq_base const& a, vector_mpq_base const& b) {
		if (a.d != b.d) throw std::runtime_error(
					"Cannot add vectors of unequal size");

		vector_mpq t(a);
		vector_mpq_base::add(t.v, b.v, t.d);
		return t;
	}

	vector_mpq operator- (vector_mpq_base const& a, vector_mpq_base const& b) {
		if (a.d != b.d) throw std::runtime_error(
					"Cannot subtract vectors of unequal size");

		vector_mpq t(a);
		vector_mpq_base::sub(t.v, b.v, t.d);
		return t;
	}

	vector_mpq operator- (vector_mpq_base const& v) {
		vector_mpq t(v.d);
		for (ind i = 0; i < v.d; ++i) {
			//t[i] = -v[i]
			mpq_neg(t.v[i].get_mpq_t(), v.v[i].get_mpq_t());
		}
		return t;
	}
	
	vector_mpq operator* (vector_mpq_base const& v, mpq_class c) {
		vector_mpq t(v);
		vector_mpq_base::mul(t.v, c, t.d);
		return t;
	}
	
	vector_mpq operator* (mpq_class c, vector_mpq_base const& v) {
		vector_mpq t(v);
		vector_mpq_base::mul(t.v, c, t.d);
		return t;
	}
	
	void vector_mpq_base::add(mpq_class* a, mpq_class const* b, ind d)
		{ for (ind i = 0; i < d; ++i) a[i] += b[i]; }

	void vector_mpq_base::sub(mpq_class* a, mpq_class const* b, ind d)
		{ for (ind i = 0; i < d; ++i) a[i] -= b[i]; }

	void vector_mpq_base::mul(mpq_class* v, mpq_class c, ind d)
		{ for (ind i = 0; i < d; ++i) v[i] *= c; }

	mpq_class inner_prod(vector_mpq_base const& a, vector_mpq_base const& b) {
		
		if (a.d != b.d) throw std::runtime_error(
			"Cannot take inner product of vectors of unequal size");
		
		mpq_class r, t;
		for (ind i = 0; i < a.d; ++i) {
			t = a.v[i] * b.v[i];
			r += t;
		}
		
		return r;		
	}
	
	vector_mpz vector_mpq_base::num() {
		vector_mpz r(d);
		for (ind i = 0; i < d; ++i) {
			mpz_set( r.v[i], v[i].get_num_mpz_t() ); /* r[i] = v[i].num() */
		}
		return r;
	}
	
	vector_mpz vector_mpq_base::den() {
		vector_mpz r(d);
		for (ind i = 0; i < d; ++i) {
			mpz_set( r.v[i], v[i].get_den_mpz_t() ); /* r[i] = v[i].den() */
		}
		return r;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t vector_mpq_hash::operator() (vector_mpq_base const& v) const {
		std::size_t seed = 0UL;
		for (ind i = 0; i < v.d; ++i) {
			/* combine low-order bits of numerator and denominator into 
				* hash */
			boost::hash_combine(seed, v.v[i].get_num().get_si() );
			boost::hash_combine(seed, v.v[i].get_den().get_si() );
		}
		return seed;
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	//
	// VECTOR_MPQ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpq::vector_mpq ( ind d ) : vector_mpq_base(new mpq_class[d], d) {}
	
	vector_mpq::vector_mpq ( vector_mpq_base const& that ) 
			: vector_mpq_base(new mpq_class[that.d], that.d) {
		for (ind i = 0; i < d; i++) v[i] = that.v[i]; 
	}
	
	vector_mpq::vector_mpq ( vector_mpq const& that ) 
			: vector_mpq_base(new mpq_class[that.d], that.d) {
		for (ind i = 0; i < d; i++) v[i] = that.v[i]; 
	}
	
	vector_mpq::vector_mpq ( matrix_row_mpq const& that ) 
			: vector_mpq_base(new mpq_class[that.d], that.d) {
		for (ind i = 0; i < d; i++) v[i] = that.v[i]; 
	}
	
	vector_mpq::vector_mpq ( vector_mpz const& that ) 
			: vector_mpq_base(new mpq_class[that.d], that.d) {
		for (ind i = 0; i < d; i++) {
			/* v[i] = that.v[i]; */
			mpq_set_z(v[i].get_mpq_t(), that.v[i]);
		}
	}
	
	vector_mpq::vector_mpq ( vector_mpz const& nums, mpz_class den ) 
			: vector_mpq_base(new mpq_class[nums.d], nums.d) {
		for (ind i = 0; i < d; i++) {
			/* v[i] = nums[i] / den; */
			mpz_set(v[i].get_num_mpz_t(), nums.v[i]);
			mpz_set(v[i].get_den_mpz_t(), den.get_mpz_t());
			v[i].canonicalize();
		}
	}
	
	vector_mpq::~vector_mpq() { delete[] v; }
	
	vector_mpq& vector_mpq::operator= ( vector_mpq_base const& that ) {
		
		if (v != that.v) {
			if ( d != that.d ) {
				delete[] v;
				d = that.d;
				v = new mpq_class[d];
			}
			
			for (ind i = 0; i < d; i++) v[i] = that.v[i];
		}
		
		return *this;
	}
	
	vector_mpq& vector_mpq::operator= ( vector_mpq const& that ) {
		
		if (v != that.v) {
			if ( d != that.d ) {
				delete[] v;
				d = that.d;
				v = new mpq_class[d];
			}
			
			for (ind i = 0; i < d; i++) v[i] = that.v[i];
		}
		
		return *this;
	}
	
	vector_mpq& vector_mpq::operator= ( matrix_row_mpq const& that ) {
		
		if (v != that.v) {
			if ( d != that.d ) {
				delete[] v;
				d = that.d;
				v = new mpq_class[d];
			}
			
			for (ind i = 0; i < d; i++) v[i] = that.v[i];
		}
		
		return *this;
	}
	
	vector_mpq& vector_mpq::operator= ( vector_mpz const& that ) {
		
		if ( d != that.d ) {
			delete[] v;
			d = that.d;
			v = new mpq_class[d];
		}
		
		for (ind i = 0; i < d; i++) {
			/* v[i] = that.v[i]; */
			mpq_set_z(v[i].get_mpq_t(), that.v[i]);
		}
		
		return *this;
	}
	
	vector_mpq& vector_mpq::operator+= (vector_mpq_base const& that) {
		if (d != that.d) throw std::runtime_error(
					"Cannot add vectors of unequal size");

		vector_mpq_base::add(v, that.v, d);
		return *this;
	}

	vector_mpq& vector_mpq::operator-= (vector_mpq_base const& that) {
		if (d != that.d) throw std::runtime_error(
					"Cannot subtract vectors of unequal size");

		vector_mpq_base::sub(v, that.v, d);
		return *this;
	}

	vector_mpq& vector_mpq::operator*= (mpq_class c) {
		vector_mpq_base::mul(v, c, d);
		return *this;
	}

	
	////////////////////////////////////////////////////////////////////////////
	//
	// MATRIX_ROW_MPQ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	matrix_row_mpq::matrix_row_mpq ( mpq_class* v, size_type d ) 
			: vector_mpq_base(v, d) {}
	
	matrix_row_mpq& matrix_row_mpq::operator= ( vector_mpq_base const& that ) {
		
		if (v != that.v) for (ind i = 0; i < d; i++) v[i] = that.v[i];
		
		return *this;
	}
	
	matrix_row_mpq& matrix_row_mpq::operator= ( vector_mpq const& that ) {
		
		if (v != that.v) for (ind i = 0; i < d; i++) v[i] = that.v[i];
		
		return *this;
	}
	
	matrix_row_mpq& matrix_row_mpq::operator= ( matrix_row_mpq const& that ) {
		
		if (v != that.v) for (ind i = 0; i < d; i++) v[i] = that.v[i];
		
		return *this;
	}
	
	matrix_row_mpq& matrix_row_mpq::operator= ( vector_mpz const& that ) {
		
		for (ind i = 0; i < d; i++) {
			/* v[i] = that.v[i]; */
			mpq_set_z(v[i].get_mpq_t(), that.v[i]);
		}
		
		return *this;
	}
	
	matrix_row_mpq& matrix_row_mpq::operator+= (vector_mpq_base const& that) {
		vector_mpq_base::add(v, that.v, d);
		return *this;
	}

	matrix_row_mpq& matrix_row_mpq::operator-= (vector_mpq_base const& that) {
		vector_mpq_base::sub(v, that.v, d);
		return *this;
	}

	matrix_row_mpq& matrix_row_mpq::operator*= (mpq_class c) {
		vector_mpq_base::mul(v, c, d);
		return *this;
	}
	
	////////////////////////////////////////////////////////////////////////////
	//
	// VECTOR_MPZ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpz::vector_mpz ( ind d ) 
			: /* lrs_alloc_mp_vector() allocs one more than needed */ 
			v(lrs_alloc_mp_vector(d-1)), d(d) {}
	
	vector_mpz::vector_mpz ( vector_mpz const& that ) 
			: /* lrs_alloc_mp_vector() allocs one more than needed */ 
			v(lrs_alloc_mp_vector(that.d-1)), d(that.d) {
		for (ind i = 0; i < d; i++) copy(v[i], that.v[i]);
	}
	
	vector_mpz::~vector_mpz() {
		/* lrs_clear_mp_vector() clears one more than needed */
		lrs_clear_mp_vector(v, d-1);
	}
	
	vector_mpz& vector_mpz::operator= ( vector_mpz const& that ) {
		if (v != that.v) {
			lrs_clear_mp_vector(v, d);
			
			d = that.d;
			v = lrs_alloc_mp_vector(d);
			for (ind i = 0; i < d; i++) copy(v[i], that.v[i]);
		}
		return *this;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Element access operators
	////////////////////////////////////////////////////////////////////////////
	
	val_t* vector_mpz::begin() { return v; }
	
	val_t const* vector_mpz::begin() const { return v; }
	
	val_t* vector_mpz::end() { return v + d; }
	
	val_t const* vector_mpz::end() const { return v + d; }
	
	val_t& vector_mpz::operator[] ( ind i ) { return v[i]; }
	
	ind vector_mpz::size() const { return d; }
	
	const val_t& vector_mpz::operator[] ( ind i ) const { return v[i]; }
	
	////////////////////////////////////////////////////////////////////////////
	// Output operator
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< ( std::ostream& o, vector_mpz const& v ) {
		o << "[";
		for (ind i = 0; i < v.d; i++) {
			o << " " << mpz_class(v[i]);
		}
		o << " ]";
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Comparison operators
	////////////////////////////////////////////////////////////////////////////
	
	bool operator<  ( vector_mpz const& a, vector_mpz const& b ) 
		{ return compare(a,b) < 0; }
	bool operator== ( vector_mpz const& a, vector_mpz const& b ) 
		{ return compare(a,b) == 0; }
	bool operator>  ( vector_mpz const& a, vector_mpz const& b ) 
		{ return compare(a,b) > 0; }
	bool operator<= ( vector_mpz const& a, vector_mpz const& b ) 
		{ return compare(a,b) <= 0; }
	bool operator!= ( vector_mpz const& a, vector_mpz const& b ) 
		{ return compare(a,b) != 0; }
	bool operator>= ( vector_mpz const& a, vector_mpz const& b ) 
		{ return compare(a,b) >= 0; }
	
	int compare( vector_mpz const& a, vector_mpz const& b ) {
		mpz_class t; int s = 0;
		
		for (ind i = 0; i < a.d && i < b.d; ++i) {
			mpz_sub(t.get_mpz_t(), a.v[i], b.v[i]);	// t = a[i] - b[i];
			s = sgn(t);
			if (s != 0) return s;
		}
		// if it reaches here, the two are lexicographically equal up to the 
		// end of the shorter vector (and s == 0)
		
		if (a.d < b.d) s = -1; else if (a.d > b.d) s = 1; /* else s = 0; */
		
		return s;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Mathematical operators
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpq vector_mpz::rationalization() const {
		if ( zero(v[0]) ) {
			/* return mpq_vector equivalent to this one */
			return vector_mpq(*this);
		} else {
			vector_mpq norm(this->d);
			for (ind i = 0; i < d; i++) {
				/* norm[i] = v[i] / v[0]; */
				mpz_set(norm.v[i].get_num_mpz_t(), v[i]);
				mpz_set(norm.v[i].get_den_mpz_t(), v[0]);
				norm.v[i].canonicalize();
			}
			return norm;
		}
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t vector_mpz_hash::operator() (vector_mpz const& v) const {
		std::size_t seed = 0UL;
		mpz_class one(1);
		for (ind i = 0; i < v.d; ++i) {
			/* combine low-order bits of value into hash */
			boost::hash_combine(seed, mpz_get_si(v.v[i]) );
			/* combine low-order bits of denominator 1 into hash (for 
			 * vector_mpq_hash compatibility) */
			boost::hash_combine(seed, one.get_si() );
		}
		return seed;
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	//
	// MATRIX_MPQ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::matrix_mpq( ind n, ind d ) 
			: m(new mpq_class[n*d]), n(n), d(d) {}
	
	matrix_mpq::matrix_mpq( matrix_mpq const& that ) 
			: m(new mpq_class[that.n*that.d]), n(that.n), d(that.d) {
		for (ind i = 0; i < n*d; i++) m[i] = that.m[i];
	}
	
	matrix_mpq::~matrix_mpq() {
		delete[] m;
	}
	
	matrix_mpq& matrix_mpq::operator= ( matrix_mpq const& that ) {
		
		if (m != that.m) {
			if ( n*d != that.n*that.d ) {
				delete[] m;
				n = that.n;
				d = that.d;
				m = new mpq_class[n*d];
			}
			
			for (ind i = 0; i < n*d; i++) m[i] = that.m[i];
		}
		
		return *this;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Element access operators
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::iterator matrix_mpq::begin() 
		{ return matrix_mpq::iterator(m, d); }
	
	matrix_mpq::const_iterator matrix_mpq::begin() const
		{ return matrix_mpq::const_iterator(m, d); }
	
	matrix_mpq::iterator matrix_mpq::end()
		{ return matrix_mpq::iterator(m+(n*d), d); }
	
	matrix_mpq::const_iterator matrix_mpq::end() const
		{ return matrix_mpq::const_iterator(m+(n*d), d); }
	
	ind matrix_mpq::size() const { return n; }
	
	ind matrix_mpq::dim() const { return d; }
	
	matrix_row_mpq matrix_mpq::operator[] ( ind i )
		{ return matrix_row_mpq(m+(i*d), d); }
	
	matrix_row_mpq const matrix_mpq::operator[] ( ind i ) const
		{ return matrix_row_mpq(m+(i*d), d); }
	
	matrix_row_mpq matrix_mpq::row( ind i )
		{ return matrix_row_mpq(m+(i*d), d); }
	
	matrix_row_mpq const matrix_mpq::row( ind i ) const 
		{ return matrix_row_mpq(m+(i*d), d); }
	
	mpq_class& matrix_mpq::elem( ind i, ind j ) { return m[i*d+j]; }
	
	mpq_class const& matrix_mpq::elem( ind i, ind j ) const { return m[i*d+j]; }

	void matrix_mpq::swap_rows(ind i, ind j) {
		for (ind k = 0; k < d; ++k) {
			mpq_class t = m[i*d+k]; m[i*d+k] = m[j*d+k]; m[j*d+k] = t;
		}
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Output operator
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< ( std::ostream& o, matrix_mpq const& m ) {
		o << "[";
		for (ind i = 0; i < m.n; i++) {
			o << " " << m[i];
		}
		o << " ]";
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Comparison operators
	////////////////////////////////////////////////////////////////////////////
	
	bool operator<  ( matrix_mpq const& a, matrix_mpq const& b )
		{ return compare(a, b) < 0; }
	
	bool operator== ( matrix_mpq const& a, matrix_mpq const& b )
		{ return compare(a, b) == 0; }
	
	bool operator>  ( matrix_mpq const& a, matrix_mpq const& b )
		{ return compare(a, b) > 0; }
	
	bool operator<= ( matrix_mpq const& a, matrix_mpq const& b )
		{ return compare(a, b) <= 0; }
	
	bool operator!= ( matrix_mpq const& a, matrix_mpq const& b )
		{ return compare(a, b) != 0; }
	
	bool operator>= ( matrix_mpq const& a, matrix_mpq const& b )
		{ return compare(a, b) >= 0; }

	int compare( matrix_mpq const& a, matrix_mpq const& b ) {
		int s = 0;
		
		for (ind i = 0; i < a.n && i < b.n; ++i) {
			s = compare(a[i], b[i]);
			if (s != 0) return s;
		}
		// if it reaches here, the two are lexicographically equal up to the 
		// end of the shorter matrix (and s == 0)
		
		if (a.n < b.n) s = -1; else if (a.n > b.n) s = 1; /* else s = 0; */
		
		return s;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Mathematical operators
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq operator* ( matrix_mpq const& a, matrix_mpq const& b ) {

		if ( a.d != b.n ) throw std::runtime_error(
						"Matrices cannot be multiplied");

		matrix_mpq c(a.n, b.d);

		if ( a.n == 0 || b.d == 0 ) return c;

		for (ind i = 0; i < a.n; ++i) for (ind j = 0; j < b.d; ++j) {
			for (ind k = 0; k < a.d; ++k) {
				c.elem(i, j) += a.elem(i, k) * b.elem(k, j);
			}
		}

		return c;
	}

	matrix_mpq operator- ( matrix_mpq const& m ) {
		matrix_mpq t(m.n, m.d);
		for (ind i = 0; i < m.n*m.d; ++i) {
			//t[i] = -m[i]
			mpq_neg(t.m[i].get_mpq_t(), m.m[i].get_mpq_t());
		}
		return t;
	}

	matrix_mpq abs(matrix_mpq const& m) {
		matrix_mpq r(m.n, m.d);
		
		for (ind i = 0; i < m.n*m.d; ++i) r.m[i] = abs(m.m[i]);
		
		return r;
	}
	
	matrix_mpq trans(matrix_mpq const& m) {
		matrix_mpq t(m.d, m.n);
		for (ind i = 0; i < m.d; ++i) for (ind j = 0; j < m.n; ++j) {
			t.elem(i, j) = m.elem(j, i);
		}
		return t;
	}

	matrix_mpq inv(matrix_mpq const& m) {

		/* inverts using Gauss-Jordan elimination */

		ind n = m.n, d = m.d;

		if ( n != d ) throw std::runtime_error(
				"Cannot invert non-square matrix");

		matrix_mpq a(m);
		matrix_mpq b(identity_mat(n));

		/* reduce down */
		for (ind i = 0; i < n; ++i) {
			if ( a.elem(i,i) != 0 ) {
				if ( a.elem(i,i) != 1 ) {
					mpq_class div;
					/* div = 1/a[i, i] */
					mpq_inv(div.get_mpq_t(), a.elem(i, i).get_mpq_t());
					b.row(i) *= div;
					a.row(i) *= div;
				}

				for (ind j = i+1; j < n; ++j) {
					if ( a.elem(j, i) != 0 ) {
						b.row(j) -= (a.elem(j, i) * b.row(i));
						a.row(j) -= (a.elem(j, i) * a.row(i));
						if ( is_zero(a.row(j)) ) {
							throw noninvertable_matrix_error(j);
						}
					}
				}
			} else {
				ind r;
				for (r = i+1; r < n; ++r) {
					if ( a.elem(r, i) != 0 ) break;
				}
				if ( r == n ) throw noninvertable_matrix_error(r);

				b.swap_rows(i, r);
				a.swap_rows(i, r);
				--i;
			}
		}

		/* reduce up */
		for (ind i = n-1; i >= 0; --i) for (ind j = i-1; j >= 0; --j) {
			if ( a.elem(j, i) != 0 ) {
				b.row(j) -= (a.elem(j, i) * b.row(i));
				a.row(j) -= (a.elem(j, i) * a.row(i));
				if ( is_zero(a.row(j)) ) throw noninvertable_matrix_error(j);
			}
		}

		return b;
	}

	matrix_mpq lu_inv(matrix_mpq const& m) {

		ind n = m.n, d = m.d;

		if ( n != d ) throw std::runtime_error(
				"Cannot invert non-square matrix");

		/* LU-decomposition */
		matrix_mpq q(n, n);

		for (ind k = 0; k < n; ++k) {
			/* Compute row of U */
			for (ind j = k; j < n; ++j) {
				mpq_class sum(0);
				for (ind s = 0; s < k; ++s) {
					sum += q.elem(k, s) * q.elem(s, j);
				}
				q.elem(k, j) = m.elem(k, j) - sum;
			}
			/* Compute column of L */
			for (ind i = k+1; i < n; ++i) {
				mpq_class sum(0);
				for (ind s = 0; s < k; ++s) {
					sum += q.elem(i, s) * q.elem(s, k);
				}
				q.elem(i, k) = (m.elem(i, k) - sum) / q.elem(k, k);
			}
		}

		/* Inverse computation */
		matrix_mpq r(n, n);

		for (ind k = 0; k < n; ++k) {
			/* Forward-solve Ly = e_k */
			vector_mpq y(n);
			for (ind i = 0; i < n; ++i) {
				y[i] = mpq_class(i == k ? 1 : 0);
				for (ind j = 0; j < i; ++j) {
					y[i] -= q.elem(i, j) * y[j];
				}
			}
			/* Backward solve Ur_k = y */
			vector_mpq x(n);
			for (ind i = n-1; i >= 0; --i) {
				x[i] = y[i];
				for (ind j = i+1; j < n; ++j) {
					x[i] -= q.elem(i, j) * x[j];
				}
				x[i] /= q.elem(i, i);
				r.elem(i, k) = x[i];
			}
		}

std::cout << "\tlu_inv():";
for (ind i = 0; i < n; ++i) {
std::cout << "\n\t";
for (ind j = 0; j < n; ++j) std::cout << " " << r.elem(i, j);
} std::cout << std::endl;

		return r;
	}

	index_set matrix_mpq::lin_indep_rows () const {
		matrix_mpq a(*this);
		int cRow = 0;
		int pivot = 0;
		std::pair<int, int> rowSwaps[n];
		int swapCounter = n;

		while ( pivot < d && cRow < n ) {

			ind swapIn = cRow;
			while ( pivot < d ) {
				bool zeroBelow = true;
				for (swapIn = cRow; swapIn < n; ++swapIn) {
					if ( a.elem(swapIn, pivot) != 0 ) {
						zeroBelow = false; break;
					}
				}
				if ( zeroBelow ) ++pivot; else break;
			}
			if ( pivot == d ) break;

			if ( swapIn != cRow ) {
				a.swap_rows(swapIn, cRow);
				rowSwaps[--swapCounter] = std::make_pair(swapIn, cRow);
			}

			if ( a.elem(cRow, pivot) != 1 ) {
				/* div = 1/a.elem(cRow, pivot) */
				mpq_class div;
				mpq_inv(div.get_mpq_t(), a.elem(cRow, pivot).get_mpq_t());
				a.row(cRow) *= div;
			}
			for (ind row = cRow+1; row < n; ++row) {
				if ( a.elem(row, pivot) != 0 ) {
					a.row(row) -= (a.elem(row, pivot) * a.row(cRow));
				}
			}

			++cRow;
			++pivot;
		}

		for (; swapCounter < n; ++swapCounter) {
			a.swap_rows(rowSwaps[swapCounter].first,
					rowSwaps[swapCounter].second);
		}

		/** 1-indexed for conformity with LRS */
		index_set r(n+1);
		for (ind i = 0; i < n; ++i) { if ( ! is_zero(a.row(i)) ) r.set(i+1); }
		return r;
	}

	matrix_mpq matrix_mpq::restriction(index_set s) const {
		matrix_mpq r(s.count(), s.count());
		
		ind i = 0;
		for (index_set_iter iterI = lrs::begin(s); 
				iterI != lrs::end(s); ++iterI) {
			ind j = 0;
			for (index_set_iter iterJ = lrs::begin(s); 
					iterJ != lrs::end(s); ++iterJ) {
				/* correct for 1-indexed index_set_iter */
				r.elem(i,j) = elem((*iterI)-1,(*iterJ)-1);
				++j;
			}
			++i;
		}
		
		return r;
	}
	
	matrix_mpq matrix_mpq::row_restriction(index_set s) const {
		matrix_mpq r(s.count(), d);

		ind i = 0;
		for (index_set_iter iterI = lrs::begin(s);
				iterI != lrs::end(s); ++iterI) {
			r.row(i) = row((*iterI)-1);
			++i;
		}

		return r;
	}

	matrix_mpq matrix_mpq::col_restriction(index_set s) const {
		matrix_mpq r(n, s.count());

		ind i = 0;
		for (ind i = 0; i < n; ++i) {
			ind j = 0;
			for (index_set_iter iterJ = lrs::begin(s);
					iterJ != lrs::end(s); ++iterJ) {
				/* correct for 1-indexed index_set_iter */
				r.elem(i,j) = elem(i,(*iterJ)-1);
				++j;
			}
		}

		return r;
	}

	vector_mpq row_mat_mul(vector_mpq_base const& v, matrix_mpq const& m) {
		ind n = m.size(), d = m.dim();

		if ( v.size() != n ) throw std::runtime_error(
				"Cannot multiply unequally sized vector and matrix");

		vector_mpq r(d);

		for (ind i = 0; i < n; ++i) for (ind j = 0; j < d; ++j) {
			r[j] += v[i] * m.elem(i, j);
		}

		return r;
	}

	matrix_mpq identity_mat(ind n) {
		matrix_mpq r(n, n);

		for (ind i = 0; i < n; ++i) r.elem(i, i) = 1;

		return r;
	}

	////////////////////////////////////////////////////////////////////////////
	// Iterator implementation
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::iterator::iterator( mpq_class* v, ind d ) : v(v), d(d) {}
	
	matrix_mpq::const_iterator::const_iterator( mpq_class* v, ind d ) 
			: v(v), d(d) {}
	
	matrix_mpq::const_iterator::const_iterator( 
			matrix_mpq::iterator const& that ) 
			: v(that.v), d(that.d) {}
	
	matrix_row_mpq matrix_mpq::iterator::operator* () 
		{ return matrix_row_mpq(v, d); }
	
	matrix_row_mpq const matrix_mpq::const_iterator::operator* () 
		{ return matrix_row_mpq(v, d); }
	
	matrix_mpq::iterator& matrix_mpq::iterator::operator++ () 
		{ v += d; return *this; }
	
	matrix_mpq::const_iterator& matrix_mpq::const_iterator::operator++ () 
		{ v += d; return *this; }
	
	matrix_mpq::iterator matrix_mpq::iterator::operator++ (int) {
		matrix_mpq::iterator it(*this);
		v += d;
		return it;
	}
	
	matrix_mpq::const_iterator matrix_mpq::const_iterator::operator++ (int) {
		matrix_mpq::const_iterator it(*this);
		v += d;
		return it;
	}
	
	matrix_mpq::iterator& matrix_mpq::iterator::operator-- () 
		{ v -= d; return *this; }
	
	matrix_mpq::const_iterator& matrix_mpq::const_iterator::operator-- () 
		{ v -= d; return *this; }
	
	matrix_mpq::iterator matrix_mpq::iterator::operator-- (int) {
		matrix_mpq::iterator it(*this);
		v -= d;
		return it;
	}
	
	matrix_mpq::const_iterator matrix_mpq::const_iterator::operator-- (int) {
		matrix_mpq::const_iterator it(*this);
		v -= d;
		return it;
	}
	
	matrix_mpq::iterator& matrix_mpq::iterator::operator+= ( ind n ) 
		{ v += n*d; return *this; }
	
	matrix_mpq::const_iterator& matrix_mpq::const_iterator::operator+= ( ind n ) 
		{ v += n*d; return *this; }
	
	matrix_mpq::iterator matrix_mpq::iterator::operator+ ( ind n ) const 
		{ return matrix_mpq::iterator(v+(n*d), d); }
	
	matrix_mpq::const_iterator matrix_mpq::const_iterator::operator+ ( 
			ind n ) const 
		{ return matrix_mpq::const_iterator(v+(n*d), d); }
	
	matrix_mpq::iterator& matrix_mpq::iterator::operator-= ( ind n ) 
		{ v -= n*d; return *this; }
	
	matrix_mpq::const_iterator& matrix_mpq::const_iterator::operator-= ( ind n ) 
		{ v -= n*d; return *this; }
	
	matrix_mpq::iterator matrix_mpq::iterator::operator- ( ind n ) const 
		{ return matrix_mpq::iterator(v-(n*d), d); }
	
	matrix_mpq::const_iterator matrix_mpq::const_iterator::operator- ( 
			ind n ) const 
		{ return matrix_mpq::const_iterator(v-(n*d), d); }
		
	ind operator- ( matrix_mpq::iterator const& a, 
					matrix_mpq::iterator const& b ) 
		{ return (a.v - b.v) / a.d; }
	
	ind operator- ( matrix_mpq::const_iterator const& a, 
					matrix_mpq::iterator const& b ) 
		{ return (a.v - b.v) / a.d; }
	
	ind operator- ( matrix_mpq::iterator const& a, 
					matrix_mpq::const_iterator const& b ) 
		{ return (a.v - b.v) / a.d; }
	
	ind operator- ( matrix_mpq::const_iterator const& a, 
					matrix_mpq::const_iterator const& b ) 
		{ return (a.v - b.v) / a.d; }
	
	bool operator== ( matrix_mpq::iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v == b.v; }
	
	bool operator== ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v == b.v; }
	
	bool operator== ( matrix_mpq::iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v == b.v; }
	
	bool operator== ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v == b.v; }
	
	bool operator!= ( matrix_mpq::iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v != b.v; }
	
	bool operator!= ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v != b.v; }
	
	bool operator!= ( matrix_mpq::iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v != b.v; }
	
	bool operator!= ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v != b.v; }
	
	bool operator<  ( matrix_mpq::iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v < b.v; }
	
	bool operator<  ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v < b.v; }
	
	bool operator<  ( matrix_mpq::iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v < b.v; }
	
	bool operator<  ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v < b.v; }
	
	bool operator>  ( matrix_mpq::iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v > b.v; }
	
	bool operator>  ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v > b.v; }
	
	bool operator>  ( matrix_mpq::iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v > b.v; }
	
	bool operator>  ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v > b.v; }
	
	bool operator<= ( matrix_mpq::iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v <= b.v; }
	
	bool operator<= ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v <= b.v; }
	
	bool operator<= ( matrix_mpq::iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v <= b.v; }
	
	bool operator<= ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v <= b.v; }
	
	bool operator>= ( matrix_mpq::iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v >= b.v; }
	
	bool operator>= ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::iterator const& b ) 
		{ return a.v >= b.v; }
	
	bool operator>= ( matrix_mpq::iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v >= b.v; }
	
	bool operator>= ( matrix_mpq::const_iterator const& a, 
					  matrix_mpq::const_iterator const& b ) 
		{ return a.v >= b.v; }
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t matrix_mpq_hash::operator() (matrix_mpq const& m) const {
		std::size_t seed = 0UL;
		for (ind i = 0; i < m.n*m.d; ++i) {
			/* combine low-order bits of numerator and denominator into 
			 * hash */
			boost::hash_combine(seed, m.m[i].get_num().get_si() );
			boost::hash_combine(seed, m.m[i].get_den().get_si() );
		}
		return seed;
	}
	
} /* namespace lrs */
