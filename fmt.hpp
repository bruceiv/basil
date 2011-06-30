#ifndef _FMT_HPP_
#define _FMT_HPP_

#include <iterator>
#include <sstream>

#include <boost/functional.hpp>

#include "dfs.hpp"

namespace basil {
	
	/** Inserts keys of a map into a set.
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					`first` member of it's dereferenced type
	 *  @param OutIter	An output iterator
	 *  @param begin	The beginning iterator of a map
	 *  @param end		The ending iterator of the map
	 *  @param out		The output iterator
	 */
	template <typename Iter, typename OutIter>
	static void insertKeys(Iter begin, Iter end, OutIter out) {
		while (begin != end) {
			*out = (*begin).first;
			++begin; ++out;
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
	static void printSet(std::ostream& o, Iter begin, Iter end) {
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << ", ";
			o << *begin;
			++begin;
		}
		o << "}";
	}
	
	/** Prints a set.
	 *  The set will be surrounded by '{' and '}', and have elements separated 
	 *  by ", ".
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					stream output operator for its dereferenced type
	 *  @param Func		A function to be applied to the dereferenced iterator
	 *  @param o		The output stream to print to
	 *  @param begin	The beginning iterator of the set
	 *  @param end		the ending iterator of the set
	 */
	template <typename Iter, typename Func>
	static void printSet(std::ostream& o, Iter begin, Iter end, Func f) {
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << ", ";
			o << f(*begin);
			++begin;
		}
		o << "}";
	}
	
	/** Returns a pointer to a specific fmt() overload
	 *  @param Arg		the type of the argument to the desired overload of 
	 *  				fmt() (without the const&)
	 */
	template <typename Arg>
	string (*fmt())(Arg const&);

	/** Prints a representation of its cobasis (as a set of indices). */
	static string fmt(dfs::index_set const& s) {
		std::stringstream o;
		printSet(o, lrs::begin(s), lrs::end(s));
		return o.str();
	}

	/** Prints a list of cobases. */
	static string fmt(dfs::cobasis_map const& m) {
		std::stringstream o;
		std::set<dfs::index_set> s;
		/* sort the set of cobases */
		insertKeys(m.begin(), m.end(), std::inserter(s, s.begin()));
		printSet(o, s.begin(), s.end(), fmt<dfs::index_set>() );
		return o.str();
	}

	/** prints a representation of a list of permutations */
	static std::ostream& operator<< (std::ostream& o, permutation_list const& l) {
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
	static string fmt(permutation_group const& g) {
		std::stringstream o;
		o << g.S;
		return o.str();
	}

	/** prints a representation of a list of vertices */
	static string fmt(dfs::coordinates_map const& m) {
		std::stringstream o;
		std::set<dfs::coordinates> s;
		/* sort the set of coordinates */
		insertKeys(m.begin(), m.end(), std::inserter(s, s.begin()));
		printSet(o, s.begin(), s.end());
		return o.str();
	}
	
	template <typename Arg>
	static string (*fmt())(Arg const&) {
		return static_cast<string (*)(Arg const&)>(&fmt);
	}
	
} /* namespace basil */

#endif /* _FMT_HPP_ */