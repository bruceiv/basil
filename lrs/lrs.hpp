#ifndef _LRS_HPP_
#define _LRS_HPP_

/** C++ wrapper class for an LRS instance.
 *
 *  @author Aaron Moss
 */

/* Copyright: Aaron Moss, 2012, moss.aaron@unb.ca         */

/* This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <iostream>
#include <new>
#include <stdexcept>
#include <string>

#include "clrs.hpp"
#include "cobasis.hpp"
#include "matrix.hpp"

/** Namespace for C++ LRS wrapper */
namespace lrs {
	
	/** Exception thrown for unexpected circumstances in the LRS wrapper.
	 *  The what string will describe the error.
	 */
	class lrs_error : public std::runtime_error {
	public:
		lrs_error(std::string const& whatArg) : runtime_error(whatArg) {}
	};
	
	/** Differentiates between expressions of equality, and expressions of 
	 *  equations */
	enum exp_type { eq = 0L, ge = 1L };
	
	/** Options for LRS.
	 *  The default constructor initializes the data members to their default 
	 *  values, while the methods can be used to modify the options.
	 */
	struct lrs_opts {
		
		lrs_opts() : out(&std::cout), vRepresentation(false) {}
		
		lrs_opts(lrs_opts const& o) 
			: out(o.out), vRepresentation(o.vRepresentation) {}
		
		lrs_opts& operator= (lrs_opts const& o) {
			out = o.out;
			vRepresentation = o.vRepresentation;
			return *this;
		}
		
		/** Sets output file */
		lrs_opts& withOutput(std::ostream& o)
			{ out = &o; return *this; }
		
		/** Gets output file */
		std::ostream& output()
			{ return *out; }
		
		/** Sets (or unsets) the V-representation flag  */
		lrs_opts& inVRepresentation(bool opt = true)
			{ vRepresentation = opt; return *this; }
		
		
		/** output stream to print output on [standard output]. */
		std::ostream* out;
		/** Specify the input is in vertex representation rather than halfspace 
		 *  representation [false]. Equivalent to the LRS input 
		 *  "V-representation" */
		bool vRepresentation;
	};
	
	/** C++ wrapper class for the LRS library. */
	class lrs {
	public:
		/** Constructor / initializer.
		 *  @param m 			the matrix to load into LRS
		 *  @param lin			the linearity indices of this matrix
		 *  @param o			LRS options (if unsupplied, will use default)
		 *  @throw bad_alloc	if the LRS process or matrix data structures 
		 * 						cannot be properly initialized.
		 */
		lrs(matrix_mpq const& m, index_set const& lin, lrs_opts o = lrs_opts());
		
		/** destructor */
		~lrs();
		
		/** Finds all candidate entering indices for a given leaving index. 
		 *  This will not find all valid pivots, but rather the nearest 
		 *  neighbours.
		 *  @param leave		The index to leave the cobasis
		 *  @return the set of indices that may enter the cobasis
		 */
		index_set allRatio(ind leave);
		
		/** Finds all candidate entering indices for a given leaving index in 
		 *  an arrangment. This will not find all valid pivots, but rather the 
		 *  nearest neighbours - it is implemented with a modified ratio test 
		 *  which finds the minimum positive and negative ratios.
		 *  @param leave		The index to leave the cobasis
		 *  @return the set of indices that may enter the cobasis
		 */
		index_set arrangementRatio(ind leave);
		
		/** finds the index in the basis array for a given entering index from 
		 *  the original list of constraints.
		 *  @param enter		The index to enter the active dictionary
		 *  @return the basis index for enter, -1 for none such found
		 */
		ind findBas(ind enter);
		
		/** finds the index in the cobasis array for a given leaving index from 
		 *  the original list of constraints.
		 *  @param leave		The index to leave the active dictionary
		 *  @return the cobasis index for leave, -1 for none such found
		 */
		ind findCob(ind leave);
		
		/** Gets the cobasis for a given column.
		 *  @param col			the column to get the cobasis for
		 *  @return a heap-allocated cobasis pointer which should be deleted by 
		 *			the caller.
		 */
		cobasis* getCobasis(ind col);
		
		/** Gets the first basis for DFS-ing from. */
		bool getFirstBasis();
		
		/** Gets the true dimension of the polytope represented by the stored 
		 *  dictionary. */
		ind getRealDim();
		
		/** Checks if the column indexed by col contains output. Returns said 
		 *  output if found.
		 *  @param col			the column index to check
		 *  @return a heap-allocated vector pointer which should be deleted by 
		 * 			the caller, or null for not a solution vector
		 *  @throw lrs_error on column negative or greater than dictionary 
		 * 			dimension.
		 */
		vector_mpz* getSolution(ind col);
		
		/** Gets solution vector for current LP ( TODO verify ) 
		 *  @return a heap-allocated vector pointer which should be deleted by 
		 * 			the caller.
		 */
		vector_mpz* getVertex();
		
		/** Finds the lex min ratio. Finds the min index ratio -aig/ais, ais<0.
		 *  @param leave		The leaving column index
		 *  @return the entering column index (-1 for none such)
		 */
		ind lexRatio(ind leave);
		
		/** Pivots the internal dictionary from leave to enter
		 *  @param leave		The leaving column index
		 *  @param enter		The entering column index
		 */
		void pivot(ind leave, ind enter);
		
		/** Prints the current dictionary */
		void printDict();
		
		/** Sets the cobasis ( TODO elaborate )
		 *  @param cob			The cobasis to set
		 *  @throw lrs_error on various error conditions
		 */
		void setCobasis(index_set& cob);
		
	private:
		
		/** Initializes LRS's problem data structure appropriately for Basil. 
		 *  Derived from lrs_read_dat (very loosely)
		 *  @param Q		The data structure to set up (should be allocated)
		 *  @param n		The number of input rows (Q->m)
		 *  @param d		The dimension of the input rows (Q->n)
		 */
		void initDat(lrs_dat* Q, ind n, ind d);
		
		/** Initializes LRS's dictionary for a problem.
		 *  Derived from lrs_read_dic (very loosely)
		 *  @param Q		The problem data (should be initialized)
		 *  @param P		The dictionary to initialize (should be allocated)
		 *  @param mat		The matrix to read in
		 *  @param lin		The linearity rows in the matrix
		 */
		void initDic(lrs_dat* Q, lrs_dic* P, matrix_mpq const& mat, 
					 index_set const& lin);
		
		/** Prints the parameter to a string. Non-negative values will be 
		 *  padded with a single space. */
		std::string toString(val_t& x);
		
		/** Structure for holding static problem data */
		lrs_dat* Q;
		/** Structure for holding current dictionary and indices */
		lrs_dic* P;
		/** Matrix holding linearities */
		lrs_mp_matrix Lin;
		/** Options provided to this instance of LRS */
		lrs_opts o;
		/** Number of instances of the LRS class */
		static int nInstances;
	}; /* class lrs */
	
} /* namespace lrs */

#endif /* _LRS_HPP_ */
