#ifndef _LRU_CACHE_HPP_
#define _LRU_CACHE_HPP_

#include <list>
#include <memory>

#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

namespace lru {
	
	/** Implements a cache with a maximum size and LRU removal semantics.
	 *  Insertion and lookup should run in O(1) time for the average case.
	 *  @param T		the type to store in the cache
	 *  @param Hash		the hash functional to use for the cache (defaults to 
	 *  				boost::hash\<T\> - note that this default may not be 
	 *  				defined for all types.)
	 *  @param Pred		the equality predicate to use for the cache (defaults 
	 *  				to std::equal_to\<T\> )
	 *  @param Alloc	the allocator to use in the hash (defaults to 
	 *  				std::allocator\<T\> )
	 *  
	 *  @author Aaron Moss
	 */
	template<
			typename T, 
			typename Hash = boost::hash<T>,
			typename Pred = std::equal_to<T>,
			typename Alloc = std::allocator<T>
		>
	class cache {
	public:
		/** Constructs a cache of the given size.
		 *  @param size		The maximum size of the cache (if 0, will be set to 
		 * 					1 - defaults to 1, but that's a bad choice)
		 */
		cache(unsigned long maxSize = 1) : cache_(), index_(), size_(0), 
				maxSize_(maxSize == 0 ? 1 : maxSize) {}
		
		/** Inserts an object into the cache. This object will be the 
		 *  most-recently used object, regardless of whether it was already in 
		 *  the cache.
		 *  @param obj		The object to insert
		 *  @return true if the object was present in the cache, false otherwise
		 */
		bool insert(T const& obj) {
			map_iter ptr = cache_.find(obj);
			bool isPresent = ( ptr != cache_.end() );
			
			if ( isPresent ) {
				//cache hit - clear previous use pointer for this object
				index_.erase( index(ptr) );
				size_--;
			}
			
			//make sure cache doesn't get overfull
			if (size_ == maxSize_) remove_lru();
			
			//put object at tail of use queue
			index_.push_back(obj);
			cache_[obj] = --(index_.end());
			size_++;
			
			return isPresent;
		}
		
		/** Looks up an object in the cache. If the object is present, it will 
		 *  become the most-recently used object.
		 *  @param obj		The object to look up
		 *  @return true if the object is found, false otherwise
		 */
		bool lookup(T const& obj) {
			map_iter ptr = cache_.find(obj);
			
			if ( ptr == cache_.end() ) {
				//cache miss
				return false;
			} else {
				//cache hit - touch object
				index_.erase( index(ptr) );
				index_.push_back(obj);
				cache_[obj] = --(index_.end());
				return true;
			}
		}
		
		/** Removes an object from the cache.
		 *  @param obj		The object to remove
		 *  @return true if the object was present in the cache, false otherwise
		 */
		bool remove(T const& obj) {
			map_iter ptr = cache_.find(obj);
			
			if ( ptr == cache_.end() ) {
				//cache miss
				return false;
			} else {
				//cache hit - remove value
				cache_.quick_erase(ptr);
				index_.erase( index(ptr) );
				size_--;
				return true;
			}
		}
		
		/** @return the current size of the cache */
		unsigned long size() { return size_; }
		
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
			while (size_ > maxSize_) remove_lru();
		}
		
	private:
		typedef 
			typename std::list<T, Alloc> 
			index_list;
		typedef 
			typename index_list::iterator 
			val_ptr;
		typedef 
			typename boost::unordered_map<T, val_ptr, Hash, Pred, Alloc> 
			cache_map;
		typedef 
			typename cache_map::iterator 
			map_iter;
		
		/** Remove the least recently used element from the cache. */
		inline void remove_lru() {
			cache_.erase( index_.front() );
			index_.pop_front();
			size_--;
		}
		
		/** Gets the index pointer from a given map pointer.
		 *  @param mp		The map pointer to return the index from
		 *  @return the index pointer from that map's value
		 */
		inline val_ptr& index(map_iter& mp) { return mp->second; }
		
		/** Hashtable based cache implementation */
		cache_map cache_;
		/** The LRU ordering on the cache. Items are ordered from least to most 
		 *  recently used. */
		index_list index_;
		/** Current size of the cache */
		unsigned long size_;
		/** Maximum size of the cache */
		unsigned long maxSize_;
			
	}; /* class lru_cache */
	
} /* namespace lru */

#endif /* _LRU_CACHE_HPP_ */
