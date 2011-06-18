#ifndef _LRU_CACHE_HPP_
#define _LRU_CACHE_HPP_

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
						/* boost::multi_index::tag<hash_lookup>,*/
						boost::multi_index::identity<T>,
						Hash,
						Pred
					>,
					/* and with secondary index an insertion-order list */
					boost::multi_index::sequenced<
						boost::multi_index::tag<lru_list>
					>
				>
			>
			cache_map;
		
		typedef
			typename cache_map::iterator
			hash_iterator;
		
		/*
		typedef 
			typename cache_map::template index<hash_lookup>::type
			hash_view;
		
		typedef
			typename hash_view::iterator
			hash_iterator;
		*/
		
		typedef
			typename cache_map::template index<lru_list>::type
			list_view;
		
		typedef
			typename list_view::iterator
			list_iterator;
		
	public:
		
		/** Type of values stored in this cache */
		typedef T value_type;
		/** Cache iterator type. Iterates from least recently used to most 
		 *  recently used type. */
		typedef list_iterator iterator;
		
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
			std::pair<hash_iterator, bool> p = cache_.insert(obj);
			
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
			//hash_view h = cache_.template get<hash_lookup>();
			//hash_iterator p = h.find(obj);
			hash_iterator p = cache_.find(obj);
			
			if ( p == cache_.end() ) {
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
			hash_iterator p = cache_.find(obj);
			
			if ( p == cache_.end() ) {
				/* cache miss */
				return false;
			} else {
				/* cache hit - remove value */
				cache_.erase(p);
				return true;
			}
		}
		
		/** Gets the beginning iterator (points to the least recently used 
		 *  item). Iteration using this iterator will not modify the use order 
		 *  of the cache. */
		iterator begin() { return cache_.template get<lru_list>().begin(); }
		
		/** Gets the ending iterator (points just past the most recently used 
		 *  item). Iteration using this iterator will not modify the use order 
		 *  of the cache. */
		iterator end() { return cache_.template get<lru_list>().end(); }
		
		/** @return the current size of the cache */
		unsigned long size() { return cache_.size(); }
		
		/** @return the maximum size of the cache */
		unsigned long maxSize() { return maxSize_; }
		
		/** Change the cache's maximum size. If the new maximum is less than 
		 *  the current number of elements, will remove the least recently used 
		 *  element until the size is at the new maximum.
		 *  @param newSize		The new maximum size of the cache. Will be set 
		 * 						to 1 if 0.
		 */
		void resize(unsigned long newSize) {
			if (newSize == 0) newSize = 1;
			maxSize_ = newSize;
			while (cache_.size() > maxSize_) eraseLru();
		}
		
	private:
		
		/** Removes the LRU element from the cache */
		inline void eraseLru() {
			cache_.template get<lru_list>().pop_front();
		}
		
		/** Touches the element of the cache given by the given iterator, 
		 *  making it most recently used. */
		inline void touch(hash_iterator p) {
			T v = *p;
			cache_.erase(p);
			cache_.insert(v);
			/* list_view& l = cache_.template get<lru_list>();
			l.relocate(cache_.template project<lru_list>(p), l.end()); */
		}
		
		/** Hashtable based cache implementation */
		cache_map cache_;
		/** Maximum size of the cache */
		unsigned long maxSize_;
			
	}; /* class lru_cache */
	
} /* namespace lru */

#endif /* _LRU_CACHE_HPP_ */
