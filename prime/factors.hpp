#ifndef _PRIME_FACTORS_HPP_
#define _PRIME_FACTORS_HPP_

#include <vector>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

#include <gmpxx.h>

/** Implements prime factorization of multi-precision integers.
 *
 *  @author Aaron Moss
 */

namespace prime {

	typedef unsigned long uind;

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
	class factorizer {
	protected:
		typedef
			boost::unordered_map<mpz_class, factor_list const, mpz_class_hash>
			factor_cache;
	public:
		/** Default constructor. */
		factorizer() : primes() {
			primes.push_back( mpz_class(2) );
			primes.push_back( mpz_class(3) );
			primes.push_back( mpz_class(5) );
			primes.push_back( mpz_class(7) );
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
} /* namespace prime */

#endif /* _PRIME_FACTORS_HPP_ */
