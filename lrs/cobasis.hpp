#ifndef _COBASIS_H_
#define _COBASIS_H_

#include <cstdlib>
#include <functional>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <gmpxx.h>

#include "clrs.hpp"

namespace lrs {
	
	/** Represents a basis or cobasis of the dictionary. For an index set s, 
	 *  s[i] == true iff index i is included in the basis or cobasis.
	 */
	typedef boost::dynamic_bitset<> index_set;
	
	/** Iterator to represent an index set as a sequence of indices, rather 
	 *  than a set of boolean values.
	 */
	class index_set_iter {
	public:
		friend index_set_iter begin(index_set& s);
		friend index_set_iter end(index_set& s);
		
		/* Uses default copy, etc. constructors */
		
		index_set_iter& operator++ () { i = s->find_next(i); return *this; }
		void operator++ (int dummy) { ++(*this); }
		
		ind operator* () { return i; }
		
		bool operator== (index_set_iter const& o)
			{ return s == o.s && i == o.i; }
		bool operator!= (index_set_iter const& o)
			{ return s != o.s || i != o.i; }
		
	private:
		
		/** Initialize this iterator on the given set with a starting index of 
		 *  the first one bit in the set. Note that this does not preclude the 
		 *  possibility that there may be no one bits, in which case the 
		 *  iterator will point to the end of the series.
		 */
		index_set_iter(index_set* s) : s(s), i(s->find_first()) { }
		
		/** Initialize this iterator on the given set with a starting index
		 *  no less than the given index. 
		 *  @param s		The set to initialize the iterator on
		 *  @param ind		The index to start the iterator on. If this index 
		 * 					does not mark a value, the next such value will be 
		 * 					chosen.
		 */
		index_set_iter(index_set* s, ind i) : s(s), i(i) { 
			if ( i < s->size() && ! s->test(i) ) i = s->find_next(i);
		}
		
		/** set this iterator is defined on */
		index_set* s;
		/** index of last seen value */
		std::size_t i;
	};
	
	/** Gets the first iterator for the index set. */
	index_set_iter begin(index_set& s) { return index_set_iter(&s); }
	/** Gets the end iterator for the index set. */
	index_set_iter end(index_set& s) { return index_set_iter(&s, s.npos); }
	
	/** Functional to hash an index_set */
	class index_set_hash : std::unary_function<index_set, std::size_t> {
	public:
		/** Hash function for an index set. XOR's the blocks of the set 
		 *  together, and returns the result.
		 *  @param s		The set to hash
		 *  @return the hash value
		 */
		std::size_t operator()(index_set const& s) const {
			/* read the block range into the XOR iterator */
			boost::to_block_range(s, it);
			/* return the final value */
			return it.v.val;
		}
		
	private:
		/** Output iterator for index set blocks, which stores the XOR'd value 
		 *  of the blocks as it sees them. */
		class xor_iter {
		friend class index_set_hash;
		private:
			/** Proxy class that overloads assignment with XOR'ing internal 
			 *  state. */
			class xor_val {
			public:
				void operator= (index_set::block_type& x) { val ^= x; }
				
				std::size_t val;
			};
			
		public:
			xor_val& operator* () { return v; }
			xor_iter& operator++ () { return *this; }
			void operator++(int dummy) {}
			
		private:
			xor_val v;
		};
		
		xor_iter it;
	};
	
	/** Gets a pseudo-random index from an index set. For an index set s, 
	 *  s.test(pseudoRandomInd(s)) is guaranteed to return true, but no 
	 *  guarantee is placed on the randomness of the distribution. Performance 
	 *  of the algorithm will be O(size(s)), but more performant algorithms 
	 *  will be preferred to more mathematically random algorithms.
	 */
	ind pseudoRandomInd(index_set& s) {
		/* generate a random index in the range [ 0 , s.size() ) */
		std::size_t randInd = std::rand() * ( s.size() - 1 ) / RAND_MAX;
		if ( ! s.test(randInd) ) {
			/* if this is not an element of the set, find the next element */
			randInd = s.find_next(randInd);
			/* if this overflowed the set, find the first element */
			if ( randInd >= s.size() ) randInd = s.find_first();
		}
		/* note that this will be biased toward indices with large empty ranges 
		 * preceding them, and is likely still worst-case linear, but should be 
		 * acceptably distributed. */
		return randInd;
	}

	
	/** Stores a cobasis and related information.
	 */
	struct cobasis {
		
		cobasis(mpz_class& det, ind ray, index_set& cob, ind totalInc, 
				index_set& extraInc) : det(det), ray(ray), cob(cob), 
				totalInc(totalInc), extraInc(extraInc) {}
		
		mpz_class det;
		ind ray;
		index_set cob;
		ind totalInc;
		index_set extraInc;
	};
	
}

#endif /* _COBASIS_H_ */
