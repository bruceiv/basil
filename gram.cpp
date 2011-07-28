/** Computes and manipulates "gram" (normed inner product representative) 
 *  matrices for the purpose of symmetry calculations.
 * 
 *  @author Aaron Moss
 */

#include <algorithm>
#include <cstdlib>
#include <ostream>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include "gram.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/matrix.hpp"

#include <iostream>

namespace basil {
	
	gram_matrix::gram_matrix(long n) : n(n), m(new int*[n]), m_(new int[n*n]) { 
		for (long i = 0; i < n; ++i) m[i] = m_+(i*n);
		std::fill(m_, m_+(n*n), 0); /* zerofill matrix */
	}
	
	gram_matrix::gram_matrix(gram_matrix const& that) 
			: n(that.n), m(new int*[that.n]), m_(new int[that.n*that.n]) { 
		for (long i = 0; i < n; ++i) /* copy row indices */
			m[i] = m_+(that.m[i]-that.m_);
		std::copy(that.m_, that.m_+(n*n), m_); /* copy matrix */
	}
		
	gram_matrix::~gram_matrix() {
		delete[] m_;
		delete[] m;
	}
	
	gram_matrix& gram_matrix::operator= (gram_matrix const& that) {
		if (m_ != that.m_) {
			if (n != that.n) {
				n = that.n;
				delete[] m_;
				delete[] m;
				m = new int*[n];
				m_ = new int[n*n];
			}
			
			for (long i = 0; i < n; ++i) /* copy row indices */
				m[i] = m_+(that.m[i]-that.m_);
			std::copy(that.m_, that.m_+(n*n), m_); /* copy matrix */
		}
		
		return *this;
	}
	
	/** Stores a multiprecision radical fraction (that is, a number of the form 
	 *  n*sqrt(r)/d), for use in gram matrix construction. */
	struct mpr {
		
		/** Constructor; initializes the mpr to canonical form [0*sqrt(1)/1] */
		mpr() : n(0), r(1), d(1) {}
		
		/** Numerator */
		mpz_class n;
		/** Radical */
		mpz_class r;
		/** Denominator */
		mpz_class d;
	};
	
	/** Equality operator for mpr type */
	bool operator== (mpr const& a, mpr const& b) 
		{ return a.n == b.n && a.r == b.r && a.d == b.d; }
	
	/** Functional to hash an mpr. */
	class mpr_hash : public std::unary_function<mpr, std::size_t> {
	public:
		std::size_t operator() (mpr const& x) const {
			std::size_t seed = 0UL;
			/* combine low-order bits of values into hash */
			boost::hash_combine(seed, x.n.get_si());
			boost::hash_combine(seed, x.r.get_si());
			boost::hash_combine(seed, x.d.get_si());
			return seed;
		}
	}; /* class mpr_hash */
	
	gram_matrix constructGram(lrs::matrix_mpq const& m) {
		
		/* typedefs */
		typedef boost::unordered_map<mpr, int, mpr_hash> mpr_map;
		
		/* gram matrix being generated */
		gram_matrix g(m.size());
		
		/* current maximum inner product representative */
		int mRep = 0;
		/* map of unique mprs to their gram matrix representatives */
		mpr_map reps;
		/* mpr representing zero */
		mpr zero;
		/* and its unique representation */
		reps.insert(std::make_pair(zero, 0));
		
		/* calculate inner products */
		mpq_class t;
		for (long i = 0; i < m.size(); ++i) {
			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (long j = 0; j <= i; ++j) {
				
// 				if ( i == 0 && j == 0 ) std::cout << "\nGRAM CONSTRUCTION";
// 				if ( j == 0 ) std::cout << std::endl;
// 				std::cout << " (" << i << "," << j << ")";
// 				std::cout << " " << m[i] << "." << m[j];
				
				t = lrs::inner_prod(m[i], m[j]);
				
				/* inner product, representing angle between m[i] and m[j], 
				 * normalized to account for rescaling of input vectors. Note 
				 * that this view of it is normalized to be positive, while the 
				 * sign is restored for the actual matrix, ensuring that the 
				 * representatives for complementary angles are negations of 
				 * each other */
				/* TODO actually do that normalization */
				mpr ip; 
				ip.n = abs(t.get_num()); 
				ip.d = t.get_den();
				
				/* look for representative for this mpr, create new if not 
				 * found */
				int rep;
				mpr_map::iterator res = reps.find(ip);
				if ( res == reps.end() ) {
					/* representative not found, make new */
					rep = ++mRep;
					reps.insert(std::make_pair(ip, rep));
				} else {
					/* representative found, use */
					rep = res->second;
				}
				rep *= sgn(t); /* sign rep to match the inner product */
				
// 				std::cout << " n" << ip.n << " r" << ip.r << " d" << ip.d 
// 						  << " -> " << rep;
				
				/* place representative in gram matrix */
				g.m[i][j] = rep; g.m[j][i] = rep;
			}
		}
		
		return g;
	}
	
	gram_matrix& gram_matrix::abs() {
		
		for (long i = 0; i < n*n; ++i) m_[i] = std::abs(m_[i]);
		
		return *this;
	}
	
	long gram_matrix::dim() const { return n; }
	
	gram_matrix gram_matrix::restriction(lrs::index_set s) const {
		gram_matrix r(s.count());
		
		long i = 0;
		for (lrs::index_set_iter iterI = lrs::begin(s); 
				iterI != lrs::end(s); ++iterI) {
			long j = 0;
			for (lrs::index_set_iter iterJ = lrs::begin(s); 
					iterJ != lrs::end(s); ++iterJ) {
				/* correct for 1-indexed index_set_iter */
				r.m[i][j] = m[(*iterI)-1][(*iterJ)-1];
				++j;
			}
			++i;
		}
		
		return r;
	}
	
	/** Functional to lexicograpically compare two rows of a gram matrix */
	class gram_matrix_row_comparator 
			: public std::binary_function<int*, int*, bool> {
	public:
		gram_matrix_row_comparator(long n = 0) : n(n) {}
		
		bool operator() (int* a, int* b) const 
			{ return std::lexicographical_compare(a, a+n, b, b+n); }
	
	private:
		/* length of the rows to compare */
		long n;
	}; /* class gram_matrix_row_comparator */
	
	gram_matrix& gram_matrix::sort() {
		
		for (long i = 0; i < n; ++i) std::sort(m[i], m[i]+n); /* sort rows */
		
		gram_matrix_row_comparator compareRows(n);
		std::sort(m, m+n, compareRows); /* lex-sort matrix by rows */
		
		return *this;
	}
	
	bool operator== (gram_matrix const& a, gram_matrix const& b) { 
		bool isEqual = a.n == b.n;
		for (long i = 0; isEqual && i < a.n; ++i) {
			isEqual = std::equal(a.m[i], a.m[i]+a.n, b.m[i]);
		}
		return isEqual;
	}
	
	std::ostream& operator<< (std::ostream& o, gram_matrix const& g) {
		for (long i = 0; i < g.n; ++i) {
			o << "| ";
			for (long j = 0; j < g.n; ++j) o << g.m[i][j] << " ";
		}
		o << "|";
		
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t gram_matrix_hash::operator() (gram_matrix const& m) const {
		std::size_t seed = 0UL;
		for (long i = 0; i < m.n; ++i) for (long j = 0; j < m.n; ++j) {
			/* combine values into hash */
			boost::hash_combine(seed, m.m[i][j] );
		}
		return seed;
	}
	
} /* namespace basil */
