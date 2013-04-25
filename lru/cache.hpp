#ifndef _LRU_CACHE_HPP_
#define _LRU_CACHE_HPP_

/** Generic least-recently-used cache data structure.
 *
 *  @author Aaron Moss
 */

/*  Copyright: Aaron Moss, 2012, moss.aaron@unb.ca  */

/*  This file is part of Basil.

    Basil is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Basil is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Basil.  If not, see <http://www.gnu.org/licenses/>.  */

#include <utility>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/sequenced_index.hpp>

namespace lru {
	
	/** Implements a cache with a maximum size and LRU removal semantics.
	 *  Insertion and lookup should run in O(1) time for the average case.
	 *  @param T		the type to store in the cache
	 *  @param Hash		the hash functional to use for the cache (defaults to 
	 *  				boost::hash\<T\> - note that this default may not be 
	 *  				defined for all types.)
	 *  @param Pred		the equality predicate to use for the cache (defaults 
	 *  				to std::equal_to\<T\> )
	 *  
	 *  @author Aaron Moss
	 */
	template<
			typename T, 
			typename Hash = boost::hash<T>,
			typename Pred = std::equal_to<T>
		>
	class cache {
	private:
		/* struct hash_lookup {}; */
		struct lru_list {};
		
		typedef
			typename boost::multi_index::multi_index_container<
				/* a container of T */
				T,
				boost::multi_index::indexed_by<
					/* with primary index a hash table with the given hash 
					 * function and equality predicate, tagged hash_lookup */
					boost::multi_index::hashed_unique<
						/* boost::multi_index::tag<hash_lookup>, */
						boost::multi_index::identity<T>,
						Hash,
						Pred
					>,
					/* with secondary index an insertion-order list, tagged 
					 * lru_list */
					boost::multi_index::sequenced<
						boost::multi_index::tag<lru_list>
					>
				>
			>
			cache_map;
		
		
		typedef
			typename cache_map::iterator
			hash_iterator;
		
		typedef
			typename cache_map::template index<lru_list>::type::iterator
			list_iterator;
		
		typedef
			hash_iterator
			primary_iterator;
		
		/*
		typedef
			typename cache_map::template index<hash_lookup>::type::iterator
			hash_iterator;
		
		
		typedef
			typename cache_map::iterator
			list_iterator;
		
		typedef
			list_iterator
			primary_iterator;
		*/
		
	public:
		
		/** Type of values stored in this cache */
		typedef T value_type;
		/** Cache iterator type. Iterates from least recently used to most 
		 *  recently used type. */
		typedef list_iterator iterator;
		/* typedef typename cache_map::reverse_iterator iterator; */
		
		/** Constructs a cache of the given size.
		 *  @param size		The maximum size of the cache (default 0)
		 */
		cache(unsigned long maxSize = 0) : cache_(), maxSize_(maxSize) {}
		
		/** Inserts an object into the cache. This object will be the 
		 *  most-recently used object, regardless of whether it was already in 
		 *  the cache.
		 *  @param obj		The object to insert
		 *  @return true if the object was present in the cache, false otherwise
		 */
		bool insert(T const& obj) {
			std::pair<primary_iterator, bool> p = add(obj);
			
			if (p.second) {
				/* Cache miss - item was not already in cache. 
				 * Ensure that cache does not exceed maximum size */
				if ( cache_.size() > maxSize() ) eraseLru();
				
			} else {
				/* Cache hit - item was already in cache.
				 * Move to most recently used. */
				touch(p.first);
			}
			
			return ! p.second;
		}
		
		/** Looks up an object in the cache. If the object is present, it will 
		 *  become the most-recently used object.
		 *  @param obj		The object to look up
		 *  @return true if the object is found, false otherwise
		 */
		bool lookup(T const& obj) {
			hash_iterator p = find(obj);
			
			if ( invalid(p) ) {
				/* cache miss */
				return false;
			} else {
				/* cache hit - touch object */
				touch(p);
				return true;
			}
		}
		
		/** Removes an object from the cache.
		 *  @param obj		The object to remove
		 *  @return true if the object was present in the cache, false otherwise
		 */
		bool remove(T const& obj) {
			hash_iterator p = find(obj);
			
			if ( invalid(p) ) {
				/* cache miss */
				return false;
			} else {
				/* cache hit - remove value */
				erase(p);
				return true;
			}
		}
		
		/** Gets the beginning iterator (points to the least recently used 
		 *  item). Iteration using this iterator will not modify the use order 
		 *  of the cache. */
		iterator begin() const {
			return cache_.template get<lru_list>().begin();
			/* return cache_.rbegin(); */
		}
		
		/** Gets the ending iterator (points just past the most recently used 
		 *  item). Iteration using this iterator will not modify the use order 
		 *  of the cache. */
		iterator end() const {
			return cache_.template get<lru_list>().end();
			/* return cache_.rend(); */
		}
		
		/** @return the current size of the cache */
		unsigned long size() const { return cache_.size(); }
		
		/** @return the maximum size of the cache */
		unsigned long maxSize() const { return maxSize_; }
		
		/** Change the cache's maximum size. If the new maximum is less than 
		 *  the current number of elements, will remove the least recently used 
		 *  element until the size is at the new maximum.
		 *  @param newSize		The new maximum size of the cache.
		 */
		void resize(unsigned long newSize) {
			maxSize_ = newSize;
			while (cache_.size() > maxSize_) eraseLru();
		}
		
	private:
		
		/** Adds an item to the cache as most recently used */
		inline std::pair<primary_iterator, bool> add(T const& obj) {
			return cache_.insert(obj);
			/* return cache_.push_front(obj); */
		}
		
		/** Removes the given element from the cache */
		inline void erase(hash_iterator p) {
			cache_.erase(p);
			/* cache_.template get<hash_lookup>().erase(p); */
		}
		
		/** Removes the LRU element from the cache */
		inline void eraseLru() {
			cache_.template get<lru_list>().pop_front();
			/* cache_.pop_back(); */
		}
		
		/** Finds an object in the cache, returning a valid hash iterator if 
		 *  it is present */
		inline hash_iterator find(T const& obj) const {
			return cache_.find(obj);
			/* return cache_.template get<hash_lookup>().find(obj); */
		}
		
		/** Tests if the given iterator is invalid */
		inline bool invalid(hash_iterator p) const {
			return p == cache_.end();
			/* return p == cache_.template get<hash_lookup>().end(); */
		}
		
		/** Touches the element of the cache given by the given iterator, 
		 *  making it most recently used. */
		/* inline void touch(list_iterator p) {
			touch(cache_.template project<hash_lookup>(p));
		} */
		
		/** Touches the element of the cache given by the given iterator, 
		 *  making it most recently used. */
		inline void touch(hash_iterator p) {
			T v = *p;
			erase(p);
			add(v);
		}
		
		
		/** Hashtable based cache implementation */
		cache_map cache_;
		/** Maximum size of the cache */
		unsigned long maxSize_;
			
	}; /* class lru_cache */
	
} /* namespace lru */

#endif /* _LRU_CACHE_HPP_ */
