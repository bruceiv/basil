#ifndef _MATRIX_GEN_HPP_
#define _MATRIX_GEN_HPP_

#include <istream>
#include <sstream>
#include <string>

#include <gmpxx.h>

#include "basilCommon.hpp"

namespace basil {
	
	/** Implementation details, not for user code. */
	namespace detail {
		
		/** Gets a line that isn't a comment (begins with '*' or '#').
		 *  @param in		the input stream to read the line from
		 *  @param s		the string to store the line in
		 *  @return the input stream
		 */
		std::istream& getContentLine(std::istream& in, string& s) {
			do {
				std::getline(in, s);
			} while (s.empty() || s[0] == '*' || s[0] == '#');
			
			return in;
		}
		
		/** Reads a line into a supplied string buffer
		 *  @param in		the input stream to read the line from
		 *  @param s		the buffer to read into
		 *  @return the original input stream 
		 */
		std::istream& getLine(std::istream& in, std::istringstream& s) {
			string t;
			std::getline(in, t);	/* read input line into temporary */
			s.clear();				/* clear error flags on stream */
			s.str(t);				/* copy new string into string buffer */
			
			return in;
		}
	} /* namespace detail */
	
	/** Allocates new matrix on heap and returns it.
	 *  Expects input in the following format (to match lrs):
	 *  
	 *  [name]
	 *  [\<lrs options\>]
	 *  begin
	 *  \<n\> \<d\> rational
	 *  \< n * d whitespace-delimited data values \>
	 *  end
	 *  
	 *  where [] denotes an optional value, \<\> denotes a variable, and any 
	 *  line beginning with '#' or '*' is ignored as a comment line.
	 */
	matrix_ptr genMatrixFromStream(std::istream& in) {
		using detail::getContentLine;
		using detail::getLine;
		
		string s = "";
		
		/* ignore lines up to begin line */
		while ( s != string("begin") ) getContentLine(in, s);
		
		/* get dimension line */
		getContentLine(in, s);
		std::istringstream lin(s);
		
		ind n, d;
		lin >> n;
		lin >> d;
		
		matrix_ptr m(new matrix(n, d));

		mpq_class t;
		for (ind i = 0; i < n; i++) {
			for (ind j = 0; j < d; j++) {
				in >> t;
				(*m)[i][j] = t;
			}
		}
		
		/* ignore up to end line */
		getContentLine(in, s);
		while ( s != string("end") ) getContentLine(in, s);
		
		return m;
	}
	
} /* namespace basil */

#endif /* _MATRIX_GEN_HPP_ */