#ifndef _PARSE_HPP_
#define _PARSE_HPP_

#include <istream>
#include <vector>

#include "basil.hpp"
#include "gram.hpp"

namespace basil {
	
	/** Input representation (constraint matrix interpretation) */
	enum representation {
		/** each row defines a halfspace in R^(d-1) */
		halfspace,
		/** each row defines a vertex in R^d */
		vertex,
		/** each row defines a line in R^(d-1) */
		arrangement
	};
	
	/** gram matrix state */
	enum gram_state { 
		/** not mentioned */
		ommited,
		/** gram matrix explicitly provided */
		provided,
		/** exact auto-generation requested */
		exact,
		/** inexact auto-generation requested */
		inexact
	};
	
	/** Stores the results of an input parse. */
	struct parse_results {
		
		/** Default constructor */
		parse_results() 
				: name(), rep(halfspace), m(), g(), l(), gm(), gs(ommited), 
				preLines(), postLines() {}
		
		/** The name of the input file */
		string name;
		/** The representation the matrix is in */
		representation rep;
		
		/** Pointer to the matrix for the problem */
		matrix_ptr m;
		/** Pointer to the permutation group for the problem */
		permutation_group_ptr g;
		/** linearity indices */
		index_set_ptr l;
		/** gram matrix for constraints */
		gram_matrix_ptr gm;
		/** parsed state of gram matrix */
		gram_state gs;
		
		/** pre-matrix unparsed lines */
		std::vector<string> preLines;
		/** post-matrix unparsed lines */
		std::vector<string> postLines;
	};
	
	/** Parses the input stream into the returned results object.
	 *  
	 *  Input is expected in the following format (based on lrs):
	 *  
	 *  [name]
	 *  [(H|V|A)-representation]
	 *  [linearity \<k\> \< k linearity indices \>]
	 *  [\< other lrs options \>]
	 *  begin
	 *  \<n\> \<d\> rational
	 *  \< n * d whitespace-delimited data values \>
	 *  end
	 *  [symmetry begin
	 *   {\<comma-delimeted cycles of whitespace delimeted elements\>}
	 *   symmetry end]
	 *  [gram auto
	 *  |gram inexact
	 *  |gram begin
	 *   \< n * n whitespace-delimited integers \>
	 *   gram end]
	 *  [\< other lrs options \>]
	 *  
	 *  @param in		The input stream to parse (will be read until its end, 
	 *  				but not closed)
	 *  @return the parse_results object encoding the information taken from 
	 *  		the stream
	 */
	parse_results parse(std::istream& in);
	
	/** Parses a matrix using the format from parse().
	 *  @param in		The input stream to parse (The "begin" line should 
	 *  				have already been consumed, and this will be consumed 
	 *  				up to the "end" line)
	 *  @return a pointer to the parsed matrix 
	 */
	matrix_ptr parseMatrix(std::istream& in);
	
	/** Parses a permutation group using the format from parse().
	 *  @param in		The input stream to parse (The "symmetry begin" line 
	 *  				should have already been consumed, and this will be 
	 *  				consumed up to the "symmetry end" line)
	 *  @param n		The number of elements the permutation group is over 
	 *  				(the number of constraints in the constraint matrix)
	 *  @return a pointer to the parsed permutation group
	 */
	permutation_group_ptr parsePermutationGroup(std::istream& in, ind n);
	
	/** Parses a gram matrix using the format from parse().
	 *  @param in		The input stream to parse (The "gram begin" line should 
	 *  				have already been consumed, and this will be consumed 
	 *  				up to the "gram end" line)
	 *  @param n		Builds an n*n gram matrix
	 *  @return a pointer to the parsed gram matrix
	 */
	gram_matrix_ptr parseGram(std::istream& in, ind n);
	
} /* namespace basil */

#endif /* _PARSE_HPP_ */
