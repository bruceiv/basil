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

#include "basil.hpp"
#include "gram.hpp"
#include "metric.hpp"

namespace basil {
	
	gram_matrix::gram_matrix(uind n, uind k_) : n(n), k_(k_), m(new int*[n]), 
			m_(new int[n*n]) { 
		for (uind i = 0; i < n; ++i) m[i] = m_+(i*n);
		std::fill(m_, m_+(n*n), 0); /* zerofill matrix */
	}
	
	gram_matrix::gram_matrix(gram_matrix const& that) 
			: n(that.n), k_(that.k_), m(new int*[that.n]), 
			m_(new int[that.n*that.n]) { 
		for (uind i = 0; i < n; ++i) /* copy row indices */
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
			
			for (uind i = 0; i < n; ++i) /* copy row indices */
				m[i] = m_+(that.m[i]-that.m_);
			std::copy(that.m_, that.m_+(n*n), m_); /* copy matrix */
			k_ = that.k_;
		}
		
		return *this;
	}
	
	int& gram_matrix::at(uind i, uind j) { return m[i][j]; }
	
	int gram_matrix::at(uind i, uind j) const { return m[i][j]; }
	
	uind gram_matrix::dim() const { return n; }
	
	gram_matrix gram_matrix::restriction(index_set s) const {
		gram_matrix r(s.count(), k_);
		
		uind i = 0;
		for (lrs::index_set_iter iterI = lrs::begin(s); 
				iterI != lrs::end(s); ++iterI) {
			uind j = 0;
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
	
	gram_matrix gram_matrix::abs() const {
		gram_matrix a(n);
		for (uind i = 0; i < n; ++i) for (uind j = 0; j < n; ++j) 
			a.m[i][j] = std::abs(m[i][j]);
		return a;
	}
	
	gram_matrix gram_matrix::doubled() const {
		gram_matrix d(2*n);
		for (uind i = 0; i < n; ++i) for (uind j = 0; j < n; ++j) {
			uind iPos = 2*i, jPos = 2*j; uind iNeg = iPos+1, jNeg = jPos+1;
			int xPos = m[i][j]; int xNeg = -xPos;
			d.m[iPos][jPos] = d.m[iNeg][jNeg] = xPos;
			d.m[iPos][jNeg] = d.m[iNeg][jPos] = xNeg;
		}
		return d;
	}
	
	gram_matrix gram_matrix::permlibCanon() const {
		
		typedef boost::unordered_map<int,int> rep_map;
		
		gram_matrix c(n);
		rep_map reps;
		
		for (uind i = 0; i < n; ++i) for (uind j = 0; j < n; ++j) {
			int val = m[i][j];
			rep_map::iterator iter = reps.find(val);
			int rep;
			if ( iter == reps.end() ) {
				rep = c.k_++;
				reps.insert(std::make_pair(val, rep));
			} else {
				rep = iter->second;
			}
			c.m[i][j] = rep;
		}
		
		return c;
	}
	
	/** Functional to lexicograpically compare two rows of a gram matrix */
	class gram_matrix_row_comparator 
			: public std::binary_function<int*, int*, bool> {
	public:
		gram_matrix_row_comparator(uind n = 0) : n(n) {}
		
		bool operator() (int* a, int* b) const 
			{ return std::lexicographical_compare(a, a+n, b, b+n); }
	
	private:
		/* length of the rows to compare */
		uind n;
	}; /* class gram_matrix_row_comparator */
	
	gram_matrix& gram_matrix::sort() {
		
		/* sort rows */
		for (uind i = 0; i < n; ++i) std::sort(m[i], m[i]+n);
		
		/* lex-sort matrix by rows */
		std::sort(m, m+n, gram_matrix_row_comparator(n));
		
		return *this;
	}
	
	bool operator== (gram_matrix const& a, gram_matrix const& b) { 
		bool isEqual = a.n == b.n;
		for (uind i = 0; isEqual && i < a.n; ++i) {
			isEqual = std::equal(a.m[i], a.m[i]+a.n, b.m[i]);
		}
		return isEqual;
	}
	
	std::ostream& operator<< (std::ostream& o, gram_matrix const& g) {
		for (uind i = 0; i < g.n; ++i) {
			o << "| ";
			for (uind j = 0; j < g.n; ++j) o << g.m[i][j] << " ";
		}
		o << "|";
		
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t gram_matrix_hash::operator() (gram_matrix const& m) const {
		std::size_t seed = 0UL;
		for (uind i = 0; i < m.n; ++i) for (uind j = 0; j < m.n; ++j) {
			/* combine values into hash */
			boost::hash_combine(seed, m.m[i][j] );
		}
		return seed;
	}
	
	/** Stores a multiprecision radical fraction (that is, a number of the form 
	 *  n*sqrt(r)/d), for use in gram matrix construction. */
	struct mpr {
		
		/** Constructor; initializes the mpr to canonical form [0*sqrt(1)/1] */
		mpr() : n(0), r(1), d(1) {}
		
		/** Constructor; initializes mpr to n_*sqrt(r_)/d_. The fraction n_/d_ 
		 *  will be reduced to lowest terms, but the radical r_ will not have 
		 *  its squares factored out. */
		mpr(mpz_class n_, mpz_class r_, mpz_class d_) : n(n_), r(r_), d(d_) {
			/* reduce the numerator and denominator to lowest terms */
			mpz_class g;
			mpz_gcd(g.get_mpz_t(), n.get_mpz_t(), d.get_mpz_t());
			n /= g; d /= g;
		}
		
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
	
	/** Inequality operator for mpr type */
	bool operator!= (mpr const& a, mpr const& b) 
		{ return !(a == b); }
	
	/** Output operator for mpr type */
	std::ostream& operator<< (std::ostream& o, mpr x) {
		o << x.n;
		if ( x.r != 1 ) o << "r" << x.r;
		if ( x.d != 1 ) o << "/" << x.d;
		
		return o;
	}
	
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
	
	/** Functional to hash an mpq_class */
	class mpq_class_hash : public std::unary_function<mpq_class, std::size_t> {
	public:
		std::size_t operator() (mpq_class const& x) const {
			std::size_t seed = 0UL;
			boost::hash_combine(seed, x.get_num().get_si());
			boost::hash_combine(seed, x.get_den().get_si());
			return seed;
		}
	}; /* class mpq_class_hash */

	/** Functional to hash an mpz_class. */
	class mpz_class_hash : public std::unary_function<mpz_class, std::size_t> {
	public:
		std::size_t operator() (mpz_class const& x) const { 
			return boost::hash_value( x.get_si() );
		}
	}; /* class mpz_class_hash */
	
	/** Postive integer representation as exponents of prime factors. */
	typedef std::vector<unsigned long> factor_list;
	
	/** Functor to compute the prime factorization of multi-precision integers. 
	 *  May maintain internal state to make subsequent factorizations faster. */
	class prime_factorizer {
	protected:
		typedef 
			boost::unordered_map<mpz_class, factor_list const, mpz_class_hash> 
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
		 *  @param x		The number to factorize (should be strictly 
		 *  				positive).
		 *  @return The prime factorization of x, as a factor_list
		 *  @throws out_of_range on x <= 0
		 */
		factor_list operator() (mpz_class x) {
			if ( x <= 0 ) {
				std::ostringstream err;
				err << "Invalid argument to prime_factorizer: " << x;
				throw std::out_of_range( err.str() );
			}
			
			factor_list l;
			
			/* look for factorization in cache first */
			factor_cache::iterator hit = cache.find(x);
			if ( hit != cache.end() ) return hit->second;
			
			/* temp that can be factored down */
			mpz_class t = x;
			
			for (uind i = 0; t > 1; ++i) {
				mpz_class p_i = 
					( i == primes.size() ) ? nextPrime() : primes[i];
				
				/* divide out all factors of p_i from t, pushing back the 
				 * number of such factors divided out into the current index in 
				 * the factor list */
				l.push_back( mpz_remove(
						t.get_mpz_t(), t.get_mpz_t(), p_i.get_mpz_t() ) 
				);
			}
			
			/* store factorization in cache */
			cache.insert(std::make_pair(x, l));
			
			return l;
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
			for (uind i = primes.size(); i < l.size(); ++i) 
				nextPrime();
			/* for each factor in the list, multiply into x, the product */
			for (uind i = 0; i < l.size(); ++i) {
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
				for (uind i = 0; i < primes.size(); ++i) {
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
	
	/** Sets rop *= op.
	 *  @param rop		One factor, where the product will be stored
	 *  @param op		The other factor
	 *  @return The product of the two
	 */
	factor_list& mult(factor_list& rop, factor_list const& op) {
		uind i;
		
		/* ensure rop is large enough to take all of op's factors */
		for (i = rop.size(); i < op.size(); ++i) rop.push_back(0);
		/* multiply by summing factors */
		for (i = 0; i < op.size(); ++i) rop[i] += op[i];
		
		return rop;
	}
	
	/** Calculates the normed inner product ip*sqrt(fi*fj)/(ni*nj). If 
	 *  ignoreSign is set, returns the absolute value of this.
	 */
	mpr norm(mpq_class ip, mpz_class ni, mpz_class nj, factor_list fi, 
			factor_list fj, prime_factorizer product) {
		
		/* if ip == 0, return default (which initializes to 0, in canonical 
		 * form) */
		if ( ip == 0 ) return mpr();
		
		/* first, multiply fj by fi */
		mult(fj, fi);
		
		/* now factor all the squares out into fi. fj, at the end should only 
		 * have 0 or 1 of each prime, and fi should have half of what was 
		 * removed  */
		
		uind k;
		/* ensure fi is big enough */
		for (k = fi.size(); k < fj.size(); ++k) fi.push_back(0);
		/* perform square root */
		for (k = 0; k < fj.size(); ++k) {
			if ( fj[k] & 0x1 /* fj[k] odd */ ) {
				fi[k] = (fj[k]-1)/2;  fj[k] = 1;
			} else /* fj[k] even */ {
				fi[k] = fj[k]/2;      fj[k] = 0;
			}
		}
		
		/* having canonicalized the radical, now generate the return value in 
		 * canonical form */
		return mpr( mpz_class( abs( ip.get_num() ) * product(fi) ),
					product( fj ),
					mpz_class( ip.get_den() * ni * nj ) );
	}
	
	gram_matrix constructEuclideanGram(matrix const& m, bool normalize) {
		
		/* typedefs */
		typedef boost::unordered_map<mpr, int, mpr_hash> mpr_map;
		
		/* size of the matrix */
		uind n = m.size();
		
		/* gram matrix being generated */
		gram_matrix g(n, 1);
		
		/* prime factorization functor */
		prime_factorizer factor;
		
		/* map of unique mprs (angles) to their gram matrix representatives */
		mpr_map reps;
		/* save representatives of zero and one as 0 and 1 */
		mpr zero; reps.insert(std::make_pair(zero, 0));
		mpr one; one.n = 1; reps.insert(std::make_pair(one, 1));
		int maxRep = 1;
		
		/* 1/||m[i]|| = sqrt(a_d*a_n)/a_n -- nums[i] = a_n, facs[i] = a_n*a_d */
		std::vector<mpz_class> nums;
		std::vector<factor_list> facs;
		
		/* calculate norm information */
		mpq_class t;
		for (uind i = 0; i < n; ++i) {
			
			if ( normalize ) {
				t = lrs::inner_prod(m[i], m[i]);
				nums.push_back(t.get_num());
				/* NOTE: this assumes here that m[i] is not a zero vector - bad 
				 * things happen otherwise */
				factor_list fn = factor( t.get_num() );
				factor_list fd = factor( t.get_den() );
				facs.push_back( mult(fn, fd) );
			}
			
			/* normed inner product of a vector with itself is 1, represented 
			 * by 1 */
			g(i,i) = 1;
		}
		
		/* calculate inner products */
//std::cout << "G:\n";
		for (uind i = 0; i < n; ++i) {
			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (uind j = 0; j < i; ++j) {
				
				t = lrs::inner_prod(m[i], m[j]);
				
				/* inner product, representing angle between m[i] and m[j],
				 * normalized to account for rescaling of input vectors. Note 
				 * that this view of it is normalized to be positive, while the 
				 * sign is restored for the actual matrix, ensuring that the 
				 * representatives for complementary angles are negations of 
				 * each other */
				
				mpr ip;
				if ( cmp(t,0) != 0 ) {
					if ( normalize ) {
						ip = norm(t, nums[i], nums[j], facs[i], facs[j], 
								  factor);
					} else {
						/* NOTE this assumes that all vectors that are 
						 * symmetric to one another have the same norm */
						ip = mpr(t.get_num(), 
								 mpz_class(1),
								 t.get_den());
					}
				}
//std::cout << " " << ip;
				
				/* look for representative for this mpr, create new if not 
				 * found */
				int rep;
				mpr_map::iterator res = reps.find(ip);
				if ( res == reps.end() ) {
					/* representative not found, make new */
					rep = ++maxRep;
					reps.insert(std::make_pair(ip, rep));
				} else {
					/* representative found, use */
					rep = res->second;
				}
				/* put proper sign on representative */
				rep *= sgn( t );
				
				/* place representative in gram matrix */
				g(i,j) = rep; g(j,i) = rep;
			}
//std::cout << "\n";
		}
//std::cout << std::endl;
		
		return g;
	}
	
	gram_matrix constructQGram(matrix const& m) {
		/* typedefs */
		typedef boost::unordered_map<mpq_class, int, mpq_class_hash> mpq_map;

		/* size of the matrix */
		ind n = m.size();
		ind d = m.dim();

/* restriction to eliminate constant terms DOESN'T WORK */
//index_set rs(d+1);
//rs.set().set(0, false).set(1, false);
//matrix rm = m.colRestriction(rs);
//std::cout << "R:\n";
//for (ind i = 0; i < rm.size(); ++i) {
//for (ind j = 0; j < rm.dim(); ++j) {
//std::cout << " " << rm.elem(i, j);
//}
//std::cout << "\n";
//}
//std::cout << std::endl;

		/* gram matrix being generated */
		gram_matrix g(n, 1);

		/* get Q-matrix inverse */
//matrix qm = rm.q_mat();
//std::cout << "Q:\n";
//for (ind i = 0; i < qm.size(); ++i) {
//for (ind j = 0; j < qm.dim(); ++j) {
//std::cout << " " << qm.elem(i,j);
//}
//std::cout << "\n";
//}
//std::cout << std::endl;
		matrix q = lu_inv(q_mat(m));
//std::cout << "inv(Q):\n";
//for (ind i = 0; i < q.size(); ++i) {
//for (ind j = 0; j < q.dim(); ++j) {
//std::cout << " " << q.elem(i,j);
//}
//std::cout << "\n";
//}
//std::cout << std::endl;

		/* get unique representatives */
		mpq_map reps;
		mpq_class zero(0); reps.insert(std::make_pair(zero, 0));
		int nextRep = 1;

//std::cout << "G:\n";
		for (ind i = 0; i < n; ++i) {
			lrs::vector_mpq w = row_mat_mul(m.row(i), q);

			for (ind j = 0; j < n; ++j) {
				int rep;
				mpq_class val = inner_prod(w, m.row(j));
//std::cout << " " << val;
				mpq_class key = abs(val);

				mpq_map::iterator res = reps.find(key);
				if ( res == reps.end() ) {
					/* representative not found, make new */
					rep = nextRep++;
					reps.insert(std::make_pair(key, rep));
				} else {
					/* representative found, use */
					rep = res->second;
				}
				rep *= sgn(val);
				g(i, j) = rep;
			}
//std::cout << "\n";
		}
//std::cout << std::endl;

		return g;
	}

} /* namespace basil */
