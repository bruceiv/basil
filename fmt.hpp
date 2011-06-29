#ifndef _FMT_HPP_
#define _FMT_HPP_

#include <sstream>

#include "dfs.hpp"

namespace basil {
	
	/** Inserts keys of a map into a set.
	 *  @param Set		A type with an insert method compatible with i->first, 
	 * 					for i of type Iter
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					`first` member of it's dereferenced type
	 *  @param s		The set to insert into
	 *  @param begin	The beginning iterator of a map
	 *  @param end		The ending iterator of the map
	 */
	template <typename Set, typename Iter>
	void insertKeys(Set s, Iter begin, Iter end) {
		while (begin != end) {
			s.insert(begin->first);
			++begin;
		}
	}
	
	/** Prints a set.
	 *  The set will be surrounded by '{' and '}', and have elements separated 
	 *  by ", ".
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					stream output operator for its dereferenced type
	 *  @param o		The output stream to print to
	 *  @param begin	The beginning iterator of the set
	 *  @param end		the ending iterator of the set
	 */
	template <typename Iter>
	void printSet(std::ostream& o, Iter begin, Iter end) {
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << ", ";
			o << *begin;
			++begin;
		}
		o << "}";
	}

	/** Prints a representation of its cobasis (as a set of indices).
	 */
	string fmt(dfs::index_set const& s) {
		std::stringstream o;
		printSet(o, lrs::begin(s), lrs::end(s));
		return o.str();
	}

	/** Prints a list of cobases. */
	string fmt(dfs::cobasis_map const& m) {
		std::stringstream o;
		std::set<dfs::index_set> s;
		insertKeys(s, m.begin(), m.end()); /* sort the set of cobases */
		printSet(o, s.begin(), s.end());
		return o.str();
	}

	/** prints a representation of a list of permutations */
	std::ostream& operator<< (std::ostream& o, permutation_list const& l) {
		bool isFirst = true;
		o << "{";
		for (permutation_list::const_iterator it = l.begin();
				it != l.end(); ++it) {
			if (isFirst) isFirst = false; else o << ", ";
			o << **it;
		}
		o << "}";
		return o;
	}

	/** prints a representation of a list of permutations */
	string fmt(permutation_group const& g) {
		std::stringstream o;
		o << g.S;
		return o.str();
	}

	/** prints a representation of a list of vertices */
	string fmt(dfs::coordinates_map const& m) {
		std::stringstream o;
		std::set<dfs::coordinates> s;
		insertKeys(s, m.begin(), m.end()); /* sort the set of coordinates */
		printSet(o, s.begin(), s.end());
		return o.str();
	}
	
} /* namespace basil */

#endif /* _FMT_HPP_ */