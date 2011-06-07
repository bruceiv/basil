#ifndef _LRU_CACHE_HPP_
#define _LRU_CACHE_HPP_

#include <list>

#include <boost/unordered_map.hpp>

namespace basil {
	
	/** Implements a cache with a maximum size and LRU removal semantics.
	 *  Insertion and lookup should run in O(1) time for the average case.
	 *  @author Aaron Moss
	 */
	template<typename T>
	class lru_cache {
	public:
		/** Constructs a cache of the given size.
		 *  @param size		The maximum size of the cache (if 0, will be set to 
		 * 					1)
		 */
		lru_cache(unsigned long maxSize) : size_(0), 
				maxSize_(maxSize == 0 ? 1 : maxSize) {}
		
		/** Inserts an object into the cache. This object will be the 
		 *  most-recently used object, regardless of whether it was already in 
		 *  the cache.
		 *  @param obj		The object to insert
		 *  @return true if the object was present in the cache, false otherwise
		 */
		bool insert(T const& obj) {
			map_ptr ptr = cache_.find(obj);
			bool isPresent = ( ptr != cache_.end() );
			
			if ( isPresent ) {
				//cache hit - clear previous use pointer for this object
				index_.erase( index(ptr) );
				size_--;
			}
			
			if (size_ == maxSize_) {
				//cache full, remove LRU element
				cache_.erase(index.front());
				index_.pop_front();
				size_--;
			}
			
			//put object at tail of use queue
			index_.push_back(obj);
			cache_[obj] = index_.back();
			size_++;
			
			return isPresent;
		}
		
		/** Looks up an object in the cache. If the object is present, it will 
		 *  become the most-recently used object.
		 *  @param obj		The object to look up
		 *  @return true if the object is found, false otherwise
		 */
		bool lookup(T const& obj) const {
			map_ptr ptr = cache_.find(obj);
			
			if ( ptr == cache_.end() ) {
				//cache miss
				return false;
			} else {
				//cache hit - touch object
				index_.erase( index(ptr) );
				index_.push_back(obj);
				cache_[obj] = index_.back();
				return true;
			}
		}
		
		/** Removes an object from the cache.
		 *  @param obj		The object to remove
		 *  @return true if the object was present in the cache, false otherwise
		 */
		bool remove(T const& obj) {
			map_ptr ptr = cache_.find(obj);
			
			if ( ptr == cache_.end() ) {
				//cache miss
				return false;
			} else {
				//cache hit - remove value
				cache_.erase(ptr);
				index_.erase( index(ptr) );
				size_--;
				return true;
			}
		}
		
		/** @return the current size of the cache */
		unsigned long size() { return size_; }
		
		/** @return the maximum size of the cache */
		unsigned long maxSize() { return maxSize_; }
		
	private:
		typedef std::list<T> index_list;
		typedef index_list::iterator val_ptr;
		typedef boost::unordered_map<T, val_ptr> cache_map;
		typedef cache_map::iterator map_ptr;
		
		/** Gets the index pointer from a given map pointer.
		 *  @param mp		The map pointer to return the index from
		 *  @return the index pointer from that map's value
		 */
		inline val_ptr& index(map_ptr& mp) { return mp->second; }
		
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
	
} /* namespace basil */

#endif /* _LRU_CACHE_HPP_ */
