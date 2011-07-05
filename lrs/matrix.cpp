/** Implements lrs::matrix class from matrix.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdlib>
#include <istream>
#include <new>
#include <ostream>
#include <stdexcept>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include "matrix.hpp"


namespace lrs {
	
	////////////////////////////////////////////////////////////////////////////
	//
	// VECTOR_MPQ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpq::vector_mpq ( ind d ) 
			: v(new mpq_t[d]), d(d), selfAlloc(true) {
		for (ind i = 0; i < d; i++) mpq_init(v[i]);
	}
	
	vector_mpq::vector_mpq ( vector_mpq const& that ) 
			: v(new mpq_t[that.d]), d(that.d), selfAlloc(true) {
		for (ind i = 0; i < d; i++) {
			mpq_init(v[i]);
			/* v[i] = that.v[i]; */
			mpq_set(v[i], that.v[i]);
		}
	}
	
	vector_mpq::vector_mpq ( vector_mpz const& that ) 
			: v(new mpq_t[that.d]), d(that.d), selfAlloc(true) {
		for (ind i = 0; i < d; i++) {
			mpq_init(v[i]);
			/* v[i] = that.v[i]; */
			mpq_set_z(v[i], that.v[i]);
		}
	}
	
	vector_mpq::vector_mpq( mpq_t* v, size_type d ) 
			: v(v), d(d), selfAlloc(false) {}
	
	vector_mpq::~vector_mpq() {
		if (selfAlloc) {
			for (ind i = 0; i < d; i++) mpq_clear(v[i]);
			delete[] v;
		}
	}
	
	vector_mpq& vector_mpq::operator= ( vector_mpq const& that ) {
		if (v != that.v) {
			for (ind i = 0; i < d; i++) mpq_clear(v[i]);
			
			if ( d != that.d && selfAlloc ) {
				delete[] v;
				d = that.d;
				v = new mpq_t[d];
			}
			
			for (ind i = 0; i < d; i++) {
				mpq_init(v[i]);
				/* v[i] = that.v[i]; */
				mpq_set(v[i], that.v[i]);
			}
		}
		return *this;
	}
	
	vector_mpq& vector_mpq::operator= ( vector_mpz const& that ) {
		for (ind i = 0; i < d; i++) mpq_clear(v[i]);
		
		if ( d != that.d && selfAlloc ) {
			delete[] v;
			d = that.d;
			v = new mpq_t[d];
		}
		
		for (ind i = 0; i < d; i++) {
			mpq_init(v[i]);
			/* v[i] = that.v[i]; */
			mpq_set_z(v[i], that.v[i]);
		}
		
		return *this;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Element access operators
	////////////////////////////////////////////////////////////////////////////
	
	mpq_t* vector_mpq::begin() { return v; }
	
	mpq_t const* vector_mpq::begin() const { return v; }
	
	mpq_t* vector_mpq::end() { return v + d; }
	
	mpq_t const* vector_mpq::end() const { return v + d; }
	
	ind vector_mpq::size() const { return d; }
	
	mpq_t& vector_mpq::operator[] ( ind i ) { return v[i]; }
	
	const mpq_t& vector_mpq::operator[] ( ind i ) const { return v[i]; }
	
	////////////////////////////////////////////////////////////////////////////
	// Output operator
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< ( std::ostream& o, vector_mpq const& v ) {
		o << "[";
		for (ind i = 0; i < v.d; i++) {
			o << " " << mpq_class(v[i]);
		}
		o << " ]";
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Comparison operators
	////////////////////////////////////////////////////////////////////////////
	
	bool operator<  ( vector_mpq const& a, vector_mpq const& b ) 
		{ return a.compare(b) < 0; }
	bool operator== ( vector_mpq const& a, vector_mpq const& b ) 
		{ return a.compare(b) == 0; }
	bool operator>  ( vector_mpq const& a, vector_mpq const& b ) 
		{ return a.compare(b) > 0; }
	bool operator<= ( vector_mpq const& a, vector_mpq const& b ) 
		{ return a.compare(b) <= 0; }
	bool operator!= ( vector_mpq const& a, vector_mpq const& b ) 
		{ return a.compare(b) != 0; }
	bool operator>= ( vector_mpq const& a, vector_mpq const& b ) 
		{ return a.compare(b) >= 0; }
	
	int vector_mpq::compare( vector_mpq const& that ) const {
		mpq_class t; int s = 0;
		
		for (ind i = 0; i < d && i < that.d; ++i) {
			mpq_sub(t.get_mpq_t(), v[i], that.v[i]);	// t = v[i] - that.v[i];
			s = sgn(t);
			if (s != 0) return s;
		}
		// if it reaches here, the two are lexicographically equal up to the 
		// end of the shorter vector (and s == 0)
		
		if (d < that.d) s = -1; else if (d > that.d) s = 1; /* else s = 0; */
		
		return s;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Mathematical operators
	////////////////////////////////////////////////////////////////////////////
	
	mpq_class inner_prod( vector_mpq const& a, vector_mpq const& b ) {
		
		if (a.d != b.d) throw std::runtime_error(
			"Cannot take inner product of vectors of unequal size");
		
		mpq_class r, t;
		for (ind i = 0; i < a.d; i++) {
			mpq_mul(t.get_mpq_t(), a.v[i], b.v[i]); /* t = a[i] * b[i] */
			r += t;
		}
		
		return r;		
	}
	
	vector_mpz vector_mpq::num() {
		vector_mpz r(d);
		for (ind i = 0; i < d; ++i) {
			mpz_set( r.v[i], mpq_numref(v[i]) ); /* r[i] = v[i].num() */
		}
		return r;
	}
	
	vector_mpz vector_mpq::den() {
		vector_mpz r(d);
		for (ind i = 0; i < d; ++i) {
			mpz_set( r.v[i], mpq_denref(v[i]) ); /* r[i] = v[i].den() */
		}
		return r;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t vector_mpq_hash::operator() (vector_mpq const& v) const {
		std::size_t seed = 0UL;
		vector_mpq& v_n = const_cast<vector_mpq&>(v);
		for (vector_mpq::iterator i = v_n.begin(); i != v_n.end(); ++i) {
			/* combine low-order bits of numerator and denominator into 
				* hash */
			boost::hash_combine(seed, mpz_get_si(mpq_numref(*i)) );
			boost::hash_combine(seed, mpz_get_si(mpq_denref(*i)) );
		}
		return seed;
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
		{ return a.compare(b) < 0; }
	bool operator== ( vector_mpz const& a, vector_mpz const& b ) 
		{ return a.compare(b) == 0; }
	bool operator>  ( vector_mpz const& a, vector_mpz const& b ) 
		{ return a.compare(b) > 0; }
	bool operator<= ( vector_mpz const& a, vector_mpz const& b ) 
		{ return a.compare(b) <= 0; }
	bool operator!= ( vector_mpz const& a, vector_mpz const& b ) 
		{ return a.compare(b) != 0; }
	bool operator>= ( vector_mpz const& a, vector_mpz const& b ) 
		{ return a.compare(b) >= 0; }
	
	int vector_mpz::compare( vector_mpz const& that ) const {
		mpz_class t; int s = 0;
		
		for (ind i = 0; i < d && i < that.d; ++i) {
			mpz_sub(t.get_mpz_t(), v[i], that.v[i]);	// t = v[i] - that.v[i];
			s = sgn(t);
			if (s != 0) return s;
		}
		// if it reaches here, the two are lexicographically equal up to the 
		// end of the shorter vector (and s == 0)
		
		if (d < that.d) s = -1; else if (d > that.d) s = 1; /* else s = 0; */
		
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
				mpz_set(mpq_numref(norm[i]), v[i]);
				mpz_set(mpq_denref(norm[i]), v[0]);
				mpq_canonicalize(norm[i]);
			}
			return norm;
		}
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t vector_mpz_hash::operator() (vector_mpz const& v) const {
		std::size_t seed = 0UL;
		vector_mpz& v_n = const_cast<vector_mpz&>(v);
		for (vector_mpz::iterator i = v_n.begin(); i != v_n.end(); ++i) {
			/* combine low-order bits of value into hash */
			boost::hash_combine(seed, mptoi(*i) );
		}
		return seed;
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	//
	// MATRIX_MPQ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::matrix_mpq( ind n, ind d ) 
			: m(new mpq_t[n*d]), n(n), d(d) {
		for (ind i = 0; i < n*d; i++) mpq_init(m[i]);
	}
	
	matrix_mpq::matrix_mpq( matrix_mpq const& that ) 
			: m(new mpq_t[that.n*that.d]), n(that.n), d(that.d) {
		for (ind i = 0; i < n*d; i++) {
			mpq_init(m[i]);
			/* m[i] = that.m[i]; */
			mpq_set(m[i], that.m[i]);
		}
	}
	
	matrix_mpq::~matrix_mpq() {
		for (ind i = 0; i < n*d; i++) mpq_clear(m[i]);
		delete[] m;
	}
	
	matrix_mpq& matrix_mpq::operator= ( matrix_mpq const& that ) {
		if (m != that.m) {
			for (ind i = 0; i < n*d; i++) mpq_clear(m[i]);
			
			if ( n*d != that.n*that.d ) {
				delete[] m;
				n = that.n;
				d = that.d;
				m = new mpq_t[n*d];
			}
			
			for (ind i = 0; i < n*d; i++) {
				mpq_init(m[i]);
				/* m[i] = that.m[i]; */
				mpq_set(m[i], that.m[i]);
			}
		}
		return *this;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Element access operators
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::iterator matrix_mpq::begin() 
		{ return matrix_mpq::iterator(m, 0, d); }
	
	matrix_mpq::const_iterator matrix_mpq::begin() const
		{ return matrix_mpq::const_iterator(m, 0, d); }
	
	matrix_mpq::iterator matrix_mpq::end()
		{ return matrix_mpq::iterator(m, n, d); }
	
	matrix_mpq::const_iterator matrix_mpq::end() const
		{ return matrix_mpq::const_iterator(m, n, d); }
	
	ind matrix_mpq::size() const { return n; }
	
	ind matrix_mpq::dim() const { return d; }
	
	vector_mpq matrix_mpq::operator[] ( ind i )
		{ return vector_mpq(m+(i*d), d); }
	
	vector_mpq const matrix_mpq::operator[] ( ind i ) const
		{ return vector_mpq(m+(i*d), d); }
	
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
		{ return a.compare(b) < 0; }
	
	bool operator== ( matrix_mpq const& a, matrix_mpq const& b )
		{ return a.compare(b) == 0; }
	
	bool operator>  ( matrix_mpq const& a, matrix_mpq const& b )
		{ return a.compare(b) > 0; }
	
	bool operator<= ( matrix_mpq const& a, matrix_mpq const& b )
		{ return a.compare(b) <= 0; }
	
	bool operator!= ( matrix_mpq const& a, matrix_mpq const& b )
		{ return a.compare(b) != 0; }
	
	bool operator>= ( matrix_mpq const& a, matrix_mpq const& b )
		{ return a.compare(b) >= 0; }

	int matrix_mpq::compare( matrix_mpq const& that ) const {
		int s = 0;
		
		for (ind i = 0; i < n && i < that.n; ++i) {
			s = operator[](i).compare(that[i]);
			if (s != 0) return s;
		}
		// if it reaches here, the two are lexicographically equal up to the 
		// end of the shorter matrix (and s == 0)
		
		if (n < that.n) s = -1; else if (n > that.n) s = 1; /* else s = 0; */
		
		return s;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Mathematical operators
	////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////
	// Iterator implementation
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::iterator::iterator( mpq_t* m, ind i, ind d ) 
			: m(m), i(i), d(d) {}
	
	matrix_mpq::const_iterator::const_iterator( mpq_t* m, ind i, ind d ) 
			: m(m), i(i), d(d) {}
	
	matrix_mpq::const_iterator::const_iterator( matrix_mpq::iterator& that ) 
			: m(that.m), i(that.i), d(that.d) {}
	
	vector_mpq matrix_mpq::iterator::dereference() const
		{ return vector_mpq(m+(i*d), d); }
	
	vector_mpq const matrix_mpq::const_iterator::dereference() const
		{ return vector_mpq(m+(i*d), d); }
	
	bool matrix_mpq::iterator::equal( matrix_mpq::iterator& that ) const
		{ return m == that.m && i == that.i; }
	
	bool matrix_mpq::iterator::equal( matrix_mpq::const_iterator& that ) const
		{ return m == that.m && i == that.i; }
	
	bool matrix_mpq::const_iterator::equal( matrix_mpq::iterator& that ) const
		{ return m == that.m && i == that.i; }
	
	bool matrix_mpq::const_iterator::equal( 
			matrix_mpq::const_iterator& that ) const
		{ return m == that.m && i == that.i; }
	
	void matrix_mpq::iterator::increment() { ++i; }
	
	void matrix_mpq::const_iterator::increment() { ++i; }
	
	void matrix_mpq::iterator::decrement() { --i; }
	
	void matrix_mpq::const_iterator::decrement() { --i; }
	
	void matrix_mpq::iterator::advance( ind n ) { i += n; }
	
	void matrix_mpq::const_iterator::advance( ind n ) { i += n; }
	
	ind matrix_mpq::iterator::distance_to( matrix_mpq::iterator& that ) const 
		{ return that.i = i; }
	
	ind matrix_mpq::iterator::distance_to( 
			matrix_mpq::const_iterator& that ) const 
		{ return that.i = i; }
	
	ind matrix_mpq::const_iterator::distance_to( 
			matrix_mpq::iterator& that ) const 
		{ return that.i = i; }
	
	ind matrix_mpq::const_iterator::distance_to( 
			matrix_mpq::const_iterator& that ) const
		{ return that.i = i; }
	
} /* namespace lrs */
