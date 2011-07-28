/** Computes and manipulates "gram" (normed inner product representative) 
 *  matrices for the purpose of symmetry calculations.
 * 
 *  @author Aaron Moss
 */

#include <algorithm>
#include <cstdlib>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
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
	}; /* struct mpr */
	
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
	
	/** Functional to hash an mpz_class. */
	class mpz_class_hash : public std::unary_function<mpz_class, std::size_t> {
	public:
		std::size_t operator() (mpz_class const& x) const { 
			return boost::hash_value( x.get_si() );
		}
	}; /* class mpz_class_hash */
	
	/** Postive integer representation as exponents of prime factors. */
	typedef std::vector<mp_bitcnt_t> factor_list;
	/** pointer to a factor_list */
	typedef boost::shared_ptr<factor_list> factor_list_ptr;
	
	/** Functor to compute the prime factorization of multi-precision integers. 
	 *  May maintain internal state to make subsequent factorizations faster. */
	class prime_factorizer {
	protected:
		typedef 
			boost::unordered_map<mpz_class, factor_list_ptr, mpz_class_hash> 
			factor_cache;
	public:
		/** Default constructor. */
		prime_factorizer() : primes() {
			primes.push_back( mpz_class(2) );
			primes.push_back( mpz_class(3) );
			/* this way I can start at the last element of the primes vector 
			 * and add two to get the first candidate */
		}
		
		/** Computes prime factorizations.
		 *  @param x		The number to factorize (should be non-negative).
		 *  @return a pointer to the prime factorization of x, as a factor_list 
		 *  		(an empty pointer if x == 0)
		 *  @throws out_of_range on x < 0
		 */
		factor_list_ptr operator() (mpz_class x) {
			if ( x < 0 ) {
				std::ostringstream err;
				err << "Invalid argument to prime_factorizer: " << x;
				throw std::out_of_range( err.str() );
			} else if ( x == 0 ) {
				return factor_list_ptr();
			}
			
			/* look for factorization in cache first */
			factor_cache::iterator hit = cache.find(x);
			if ( hit != cache.end() ) return hit->second;
			
			factor_list_ptr p = boost::make_shared<factor_list>();
			factor_list& l = *p;
			
			for (unsigned long i = 0; x > 1; ++i) {
				mpz_class p_i = 
					( i == primes.size() ) ? nextPrime() : primes[i];
				
				/* divide out all factors of p_i from x, pushing back the 
				 * number of such factors divided out into the current index in 
				 * the factor list */
				l.push_back( mpz_remove(
						x.get_mpz_t(), x.get_mpz_t(), p_i.get_mpz_t() ) 
				);
			}
			
			/* store factorization in cache */
			cache.insert(std::make_pair(x, p));
			
			return p;
		}
		
		/** Converts a factor list to an integer
		 *  @param l		The list to multiply out
		 *  @return the product of the factors denoted in l, as a 
		 *  		multi-precision integer.
		 */
		mpz_class operator() (factor_list const& l) {
			
			mpz_class x(1);
			mpz_class t;
			
			/* ensure that primes is long enough to handle this list */
			for (unsigned long i = primes.size(); i < l.size(); ++i) 
				nextPrime();
			/* for each factor in the list, multiply into x, the product */
			for (unsigned long i = 0; i < l.size(); ++i) {
				switch( l[i] ) {
				case 0: /* do nothing */	break;
				case 1: x *= primes[i];		break;
				default: 
					/* t = primes[i]^l[i] (exponentiation, not XOR) */
					mpz_pow_ui(t.get_mpz_t(), primes[i].get_mpz_t(), l[i]);
					x *= t;
					break;
				}
			}
			
			return x;
		}
		
	protected:
		
		/** Adds a new prime to the end of the primes vector.
		 *  @return the new prime 
		 */
		mpz_class nextPrime() {
			/* use the Sieve of Eratsothenes to compute the next prime */
			mpz_class cand = primes.back();
			
			bool foundPrime = false;
			while ( ! foundPrime ) {
				cand += 2;
				for (unsigned long i = 0; i < primes.size(); ++i) {
					mpz_class& p_i = primes[i];
					
					int c = cmp(mpz_class(p_i*p_i), cand);
					/* candidate is less than p_i^2, therefore prime */
					if ( c > 0 ) { foundPrime = true; break; }
					/* candidate is p_i^2, therefore composite */
					else if ( c == 0 ) break;
					
					/* prime divides candidate, therefore composite */
					if ( cand % p_i == 0 ) break;
				}
			}
			
			primes.push_back(cand);
			return cand;
		}
		
		/** The primes used in factorization. */
		std::vector<mpz_class> primes;
		/** Cache of already-computed factorizations */
		factor_cache cache;
	}; /* class prime_factorizer */
	
	/** Calculates the normed inner product |ip|*sqrt(fi*fj)/(ni*nj) */
	mpr norm(mpq_class ip, mpz_class ni, mpz_class nj, factor_list fi, 
			factor_list fj, prime_factorizer product) {
		
		mpr x;
		/* if ip == 0, return x (which initializes to 0, in canonical form) */
		if ( ip == 0 ) return x;
		
		/* first, multiply fi by fj -- easy, just add factors */
		unsigned long k;
		/* ensure the factor lists are the same length */
		for (k = fi.size(); k < fj.size(); ++k) fi.push_back(0);
		for (k = fj.size(); k < fi.size(); ++k) fj.push_back(0);
		/* sum factors into fj */
		for (k = 0; k < fj.size(); ++k) fj[k] += fi[k];
		
		/* now factor all the squares out into fi - fj, at the end should only 
		 * have 0 or 1 of each prime, and fi should have half of what was 
		 * removed  */
		for (k = 0; k < fj.size(); ++k) {
			if ( fj[k] & 0x1 /* fj[k] odd */ ) {
				fi[k] = (fj[k]-1)/2;  fj[k] = 1;
			} else /* fj[k] even */ {
				fi[k] = fj[k]/2;      fj[k] = 0;
			}
		}
		
		/* initialize the return value */
		x.n = mpz_class( abs(ip.get_num()) * product(fi) );
		x.r = product( fj );
		x.d = mpz_class( ip.get_den() * ni * nj );
		
		/* now reduce its numerator and denominator to lowest terms */
		mpz_class g;
		mpz_gcd(g.get_mpz_t(), x.n.get_mpz_t(), x.d.get_mpz_t());
		x.n /= g; x.d /= g;
		
		return x;
	}
	
	gram_matrix constructGram(lrs::matrix_mpq const& m) {
		
		/* typedefs */
		typedef boost::unordered_map<mpr, int, mpr_hash> mpr_map;
		
		/* gram matrix being generated */
		gram_matrix g(m.size());
		
		/* prime factorization functor */
		prime_factorizer factor;
		
		/* current maximum inner product representative */
		int mRep = 1;
		/* map of unique mprs to their gram matrix representatives */
		mpr_map reps;
		
		/* preload 0 and 1 */
		mpr zero; reps.insert(std::make_pair(zero, 0));
		mpr one; one.n = 1; reps.insert(std::make_pair(one, 1));
		
		/* 1/||m[i]|| = sqrt(a_d*a_n)/a_n -- nums[i] = a_n, facs[i] = a_n*a_d */
		std::vector<mpz_class> nums;
		std::vector<factor_list_ptr> facs;
		
		/* calculate norm information */
		mpq_class t;
		for (long i = 0; i < m.size(); ++i) {
			t = lrs::inner_prod(m[i], m[i]);
			
			nums.push_back(t.get_num());
			facs.push_back( factor( mpz_class( t.get_num() * t.get_den() ) ) );
			
			g.m[i][i] = 1; /* inner product of any vector with itself is 1 */
		}
		
		/* calculate inner products */
		for (long i = 0; i < m.size(); ++i) {
			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (long j = 0; j < i; ++j) {
				
				if ( i == 0 && j == 0 ) std::cout << "\nGRAM CONSTRUCTION";
				if ( j == 0 ) std::cout << std::endl;
				std::cout << " (" << i << "," << j << ")";
				std::cout << " " << m[i] << "." << m[j];
				
				t = lrs::inner_prod(m[i], m[j]);
				
				/* inner product, representing angle between m[i] and m[j], 
				 * normalized to account for rescaling of input vectors. Note 
				 * that this view of it is normalized to be positive, while the 
				 * sign is restored for the actual matrix, ensuring that the 
				 * representatives for complementary angles are negations of 
				 * each other */
				
				mpr ip = norm(t, nums[i], nums[j], *facs[i], *facs[j], factor); 
// 				ip.n = abs(t.get_num()); 
// 				ip.d = t.get_den();
				
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
				
				std::cout << " n" << ip.n << " r" << ip.r << " d" << ip.d 
						  << " -> " << rep;
				
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
