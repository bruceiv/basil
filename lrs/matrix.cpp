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

#include "cobasis.hpp"
#include "matrix.hpp"


namespace lrs {
	
	////////////////////////////////////////////////////////////////////////////
	//
	// VECTOR_MPQ CLASS
	//
	////////////////////////////////////////////////////////////////////////////
	
	vector_mpq::vector_mpq ( ind d ) 
			: v(new mpq_class[d]), d(d), selfAlloc(true) {}
	
	vector_mpq::vector_mpq ( vector_mpq const& that ) 
			: v(new mpq_class[that.d]), d(that.d), selfAlloc(true) {
		for (ind i = 0; i < d; i++) v[i] = that.v[i]; 
	}
	
	vector_mpq::vector_mpq ( vector_mpz const& that ) 
			: v(new mpq_class[that.d]), d(that.d), selfAlloc(true) {
		for (ind i = 0; i < d; i++) {
			/* v[i] = that.v[i]; */
			mpq_set_z(v[i].get_mpq_t(), that.v[i]);
		}
	}
	
	vector_mpq::vector_mpq( mpq_class* v, size_type d ) 
			: v(v), d(d), selfAlloc(false) {}
	
	vector_mpq::~vector_mpq() {
		if (selfAlloc) delete[] v;
	}
	
	vector_mpq& vector_mpq::operator= ( vector_mpq const& that ) {
		
		if (v != that.v) {
			
			if ( d != that.d && selfAlloc ) {
				delete[] v;
				d = that.d;
				v = new mpq_class[d];
			}
			
			for (ind i = 0; i < d; i++) v[i] = that.v[i];
			
		}
		return *this;
	}
	
	vector_mpq& vector_mpq::operator= ( vector_mpz const& that ) {
		
		if ( d != that.d && selfAlloc ) {
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
	
	////////////////////////////////////////////////////////////////////////////
	// Element access operators
	////////////////////////////////////////////////////////////////////////////
	
	mpq_class* vector_mpq::begin() { return v; }
	
	mpq_class const* vector_mpq::begin() const { return v; }
	
	mpq_class* vector_mpq::end() { return v + d; }
	
	mpq_class const* vector_mpq::end() const { return v + d; }
	
	ind vector_mpq::size() const { return d; }
	
	mpq_class& vector_mpq::operator[] ( ind i ) { return v[i]; }
	
	const mpq_class& vector_mpq::operator[] ( ind i ) const { return v[i]; }
	
	////////////////////////////////////////////////////////////////////////////
	// Output operator
	////////////////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< ( std::ostream& o, vector_mpq const& v ) {
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
			t = v[i] - that.v[i];
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
			t = a.v[i] * b.v[i];
			r += t;
		}
		
		return r;		
	}
	
	vector_mpz vector_mpq::num() {
		vector_mpz r(d);
		for (ind i = 0; i < d; ++i) {
			mpz_set( r.v[i], v[i].get_num_mpz_t() ); /* r[i] = v[i].num() */
		}
		return r;
	}
	
	vector_mpz vector_mpq::den() {
		vector_mpz r(d);
		for (ind i = 0; i < d; ++i) {
			mpz_set( r.v[i], v[i].get_den_mpz_t() ); /* r[i] = v[i].den() */
		}
		return r;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t vector_mpq_hash::operator() (vector_mpq const& v) const {
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
	
	vector_mpq matrix_mpq::row( ind i )
		{ return vector_mpq(m+(i*d), d); }
	
	vector_mpq const matrix_mpq::row( ind i ) const 
		{ return vector_mpq(m+(i*d), d); }
	
	mpq_class& matrix_mpq::elem( ind i, ind j ) { return m[i*d+j]; }
	
	mpq_class const& matrix_mpq::elem( ind i, ind j ) const { return m[i*d+j]; }

	
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
			s = row(i).compare(that.row(i));
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
	
	matrix_mpq matrix_mpq::inner_prod_mat () {
		matrix_mpq p(n,n);
		
		mpq_class t;
		for (ind i = 0; i < n; ++i) {
			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (ind j = 0; j < i; ++j) {
				t = inner_prod(row(i), row(j));
				p.elem(i,j) = t; 
				p.elem(j,i) = t;
			}
			/* Handle j = i case here */
			t = inner_prod(row(i), row(i));
			p.elem(i,i) = t;
		}
		
		return p;
	}
	
	matrix_mpq matrix_mpq::restriction(index_set s) {
		matrix_mpq r(s.count(), s.count());
		
		ind i = 0;
		for (index_set_iter iterI = lrs::begin(s); 
				iterI != lrs::end(s); ++iterI) {
			ind j = 0;
			for (index_set_iter iterJ = lrs::begin(s); 
					iterJ != lrs::end(s); ++iterJ) {
				r.elem(i,j) = elem(*iterI,*iterJ);
				++j;
			}
			++i;
		}
		
		return r;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Iterator implementation
	////////////////////////////////////////////////////////////////////////////
	
	matrix_mpq::iterator::iterator( mpq_class* m, ind i, ind d ) 
			: m(m), i(i), d(d) {}
	
	matrix_mpq::const_iterator::const_iterator( mpq_class* m, ind i, ind d ) 
			: m(m), i(i), d(d) {}
	
	matrix_mpq::const_iterator::const_iterator( 
			matrix_mpq::iterator const& that ) 
			: m(that.m), i(that.i), d(that.d) {}
	
	vector_mpq matrix_mpq::iterator::dereference() const
		{ return vector_mpq(m+(i*d), d); }
	
	vector_mpq const matrix_mpq::const_iterator::dereference() const
		{ return vector_mpq(m+(i*d), d); }
	
	bool matrix_mpq::iterator::equal( matrix_mpq::iterator const& that ) const
		{ return m == that.m && i == that.i; }
	
	bool matrix_mpq::iterator::equal( 
			matrix_mpq::const_iterator const& that ) const
		{ return m == that.m && i == that.i; }
	
	bool matrix_mpq::const_iterator::equal( 
			matrix_mpq::iterator const& that ) const
		{ return m == that.m && i == that.i; }
	
	bool matrix_mpq::const_iterator::equal( 
			matrix_mpq::const_iterator const& that ) const
		{ return m == that.m && i == that.i; }
	
	void matrix_mpq::iterator::increment() { ++i; }
	
	void matrix_mpq::const_iterator::increment() { ++i; }
	
	void matrix_mpq::iterator::decrement() { --i; }
	
	void matrix_mpq::const_iterator::decrement() { --i; }
	
	void matrix_mpq::iterator::advance( ind n ) { i += n; }
	
	void matrix_mpq::const_iterator::advance( ind n ) { i += n; }
	
	ind matrix_mpq::iterator::distance_to( 
			matrix_mpq::iterator const& that ) const 
		{ return that.i - i; }
	
	ind matrix_mpq::iterator::distance_to( 
			matrix_mpq::const_iterator const& that ) const 
		{ return that.i - i; }
	
	ind matrix_mpq::const_iterator::distance_to( 
			matrix_mpq::iterator const& that ) const 
		{ return that.i - i; }
	
	ind matrix_mpq::const_iterator::distance_to( 
			matrix_mpq::const_iterator const& that ) const
		{ return that.i - i; }
	
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
