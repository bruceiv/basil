#ifndef _LRS_HPP_
#define _LRS_HPP_

#include <iostream>
#include <new>
#include <stdexcept>
#include <string>

#include "clrs.hpp"
#include "cobasis.hpp"
#include "matrix.hpp"

/** Namespace for C++ LRS wrapper */
namespace lrs {
	
	/** Exception thrown for unexpected circumstances in the DFS algorithm.
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
		 *  @param o			LRS options (if unsupplied, will use default)
		 *  @throw bad_alloc	if the LRS process or matrix data structures 
		 * 						cannot be properly initialized.
		 */
		lrs(matrix const& m, lrs_opts o = lrs_opts()) throw(std::bad_alloc);
		
		/** destructor */
		~lrs();
		
		/** finds all candidate entering indices for a given leaving index
		 *  @param leave		The index to leave the dictionary
		 *  @return the set of indices that may enter the dictionary
		 */
		index_set allRatio(ind leave);
		
		/** finds the basis for a given entering index.
		 *  @param enter		The index to enter the active dictionary
		 *  @return the basis for that index, -1 for none such found
		 */
		ind findBas(ind enter);
		
		/** finds the cobasis for a given leaving index.
		 *  @param leave		The index to leave the active dictionary
		 *  @return the cobasis for that index, -1 for none such found
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
		 */
		void initDic(lrs_dat* Q, lrs_dic* P, matrix const& mat);
		
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
	}; /* class lrs */
	
} /* namespace lrs */

#endif /* _LRS_HPP_ */
