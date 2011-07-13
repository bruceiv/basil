#ifndef _FMT_HPP_
#define _FMT_HPP_

#include <iterator>
#include <sstream>

#include <boost/functional.hpp>

#include "dfs.hpp"
#include "lrs/cobasis.hpp"

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
	 *  by delim.
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					stream output operator for its dereferenced type
	 *  @param o		The output stream to print to
	 *  @param begin	The beginning iterator of the set
	 *  @param end		the ending iterator of the set
	 *  @param delim	the set element delimiter [", "]
	 */
	template <typename Iter>
	static void printSet(std::ostream& o, Iter begin, Iter end, 
						 string delim = ", ") {
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << delim;
			o << *begin;
			++begin;
		}
		o << "}";
	}
	
	/** Prints a set.
	 *  The set will be surrounded by '{' and '}', and have elements separated 
	 *  by delim.
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					stream output operator for its value type
	 *  @param Func		A unary function taking an argument of the iterator's 
	 * 					value type, and returning an object with a standard 
	 * 					output operator defined
	 *  @param o		The output stream to print to
	 *  @param begin	The beginning iterator of the set
	 *  @param end		the ending iterator of the set
	 *  @param f		A function of type Func
	 *  @param delim	The set element delimeter [", "]
	 */
	template <typename Iter, typename Func>
	static void printSet(std::ostream& o, Iter begin, Iter end, Func f, 
						 string delim = ", ") {
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << delim;
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
		typedef
			bool (*index_set_comparator)(dfs::index_set const&, 
										 dfs::index_set const&);
		
		std::stringstream o;
		std::set<dfs::index_set, index_set_comparator> 
			s(&lrs::lexicographical_compare);
		/* sort the set of cobases */
		insertKeys(m.begin(), m.end(), std::inserter(s, s.begin()));
		printSet(o, s.begin(), s.end(), fmt<dfs::index_set>() );
		return o.str();
	}

	/** prints a representation of a list of permutations */
	static std::ostream& operator<< (std::ostream& o, 
									 permutation_list const& l) {
		printSet(o, l.begin(), l.end(), 
				 boost::mem_fun_ref(&permlib::Permutation::ptr::operator*) );
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