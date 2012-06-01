#ifndef _FMT_HPP_
#define _FMT_HPP_

/** Pretty-printing methods for a variety of Basil types.
 *
 *  @author Aaron Moss
 */

#include <functional>
#include <iterator>
#include <sstream>

#include <boost/functional.hpp>

#include "dfs_types.hpp"

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
	
	/** Integer constant for single line set printing */
	static const int single_line = -20; /* the big negative is to allow for 
										 * some increment. If you add more than 
										 * this many tabstops, formatting will 
										 * be ugly anyway */
	
	/** Gets the line delimiter for the given number of tab stops
	 *  @param tabs		The number of tab stops to place. If less than 0, no 
	 * 					newline, delimeter is a single space.
	 *  @return a newline followed by tabs tab characters */
	static string lineSpace(int tabs) {
		if (tabs < 0) return string(" ");
		string s(tabs+1, '\t'); s[0] = '\n';
		return s;
	}
	
	/** Prints a set.
	 *  The set will be surrounded by '{' and '}', and have elements separated 
	 *  by ','.
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					stream output operator for its dereferenced type
	 *  @param o		The output stream to print to
	 *  @param begin	The beginning iterator of the set
	 *  @param end		the ending iterator of the set
	 *  @param tabs		number of tabstops to indent set (elements will be 
	 * 					indented one stop further) [default is single line]
	 */
	template <typename Iter>
	static void printSet(std::ostream& o, Iter begin, Iter end, 
						 int tabs = single_line) {
		string space = lineSpace(tabs+1);
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << ",";
			o << space << *begin;
			++begin;
		}
		if (! isFirst ) o << lineSpace(tabs);
		o << "}";
	}
	
	/** Prints a set.
	 *  The set will be surrounded by '{' and '}', and have elements separated 
	 *  by ','
	 *  @param Iter		A type that has pre-increment defined, as well as a 
	 * 					stream output operator for its value type
	 *  @param Func		A unary function taking an argument of the iterator's 
	 * 					value type, and returning an object with a standard 
	 * 					output operator defined
	 *  @param o		The output stream to print to
	 *  @param begin	The beginning iterator of the set
	 *  @param end		the ending iterator of the set
	 *  @param tabs		number of tabstops to indent set (elements will be 
	 * 					indented one stop further) [default is single line]
	 *  @param f		A function of type Func
	 */
	template <typename Iter, typename Func>
	static void printSet(std::ostream& o, Iter begin, Iter end, Func f,
						 int tabs = single_line) {
		string space = lineSpace(tabs+1);
		bool isFirst = true;
		o << "{";
		while (begin != end) {
			if (isFirst) isFirst = false; else o << ",";
			o << space << f(*begin);
			++begin;
		}
		if (! isFirst ) o << lineSpace(tabs);
		o << "}";
	}
	
	/** Returns a pointer to a specific fmt() overload
	 *  @param Arg		the type of the argument to the desired overload of 
	 *  				fmt() (without the const&)
	 *  @param tabs		The number of tabstops to set the overload to [default 
	 *  				single line]
	 */
	template <typename Arg>
	boost::binder2nd<string (*)(Arg const&, int)> fmt(int tabs = single_line);

	/** Prints a representation of its cobasis (as a set of indices). */
	static string fmt(index_set const& s, int tabs = single_line) {
		std::ostringstream o;
		printSet(o, lrs::begin(s), lrs::end(s), tabs);
		return o.str();
	}

	/** Prints a list of cobases. */
	static string fmt(cobasis_map const& m, int tabs = single_line) {
		typedef
			bool (*index_set_comparator)(index_set const&, index_set const&);
		
		std::ostringstream o;
		std::set<index_set, index_set_comparator> 
				s(&lrs::lexicographical_compare);
		/* sort the set of cobases */
		insertKeys(m.begin(), m.end(), std::inserter(s, s.begin()));
		o << s.size() << lineSpace(tabs);
		printSet(o, s.begin(), s.end(), fmt<index_set>(), tabs);
		return o.str();
	}

	/** prints a representation of a list of permutations */
	static string fmt(permutation_group const& g, int tabs = single_line) {
		std::ostringstream o;
		o << g.S.size() << lineSpace(tabs);
		printSet(o, g.S.begin(), g.S.end(),
				 boost::mem_fun_ref(&permlib::Permutation::ptr::operator*), 
				 tabs);
		return o.str();
	}

	/** prints a representation of a list of vertices */
	static string fmt(coordinates_map const& m, int tabs = single_line) {
		std::ostringstream o;
		std::set<coordinates> s;
		/* sort the set of coordinates */
		insertKeys(m.begin(), m.end(), std::inserter(s, s.begin()));
		o << s.size() << lineSpace(tabs);
		printSet(o, s.begin(), s.end(), tabs);
		return o.str();
	}
	
	/** prints a matrix */
	static string fmt(matrix const& m, int tabs = single_line) {
		std::ostringstream o;
		string row_space = (tabs < 0) ? " " : (lineSpace(tabs) + "  ");
		o << "(" << m.size() << "," << m.dim() << ")" << lineSpace(tabs);
		
		bool isFirst = true;
		matrix::const_iterator iter = m.begin(), end = m.end();
		o << "[";
		while ( iter != end ) {
			if (isFirst) {
				isFirst = false; 
				o << " ";
			} else {
				o << row_space;
			}
			o << *iter;
			++iter;
		}
		o << " ]";
		
		return o.str();
	}
	
	template <typename Arg>
	boost::binder2nd<string (*)(Arg const&, int)> fmt(int tabs) {
		return boost::bind2nd(
			static_cast<string (*)(Arg const&, int)>(&fmt), 
			tabs);
	}
	
} /* namespace basil */

#endif /* _FMT_HPP_ */
