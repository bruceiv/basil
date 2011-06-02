#ifndef _LRS_HPP_
#define _LRS_HPP_

#include <new>

#include "clrs.hpp"
#include "matrix.hpp"

/** Namespace for C++ LRS wrapper */
namespace lrs {
	
	/** Differentiates between expressions of equality, and expressions of 
	 *  equations */
	enum exp_type { eq = 0L, ge = 1L };
	
	/** C++ wrapper class for the LRS library. */
	class lrs {
	public:
		/** Constructor / initializer.
		 *  @param m 			the matrix to load into LRS
		 *  @throw bad_alloc	if the LRS process or matrix data structures 
		 * 						cannot be properly initialized.
		 */
		lrs(const matrix& m) throw(std::bad_alloc);
		
		/** destructor */
		~lrs();
		
		/** Gets the first basis for DFS-ing from. */
		bool getFirstBasis();
		
		/** Prints the current dictionary */
		void printDict();
		
	private:
		
		/** Structure for holding static problem data */
		lrs_dat* Q;
		/** Structure for holding current dictionary and indices */
		lrs_dic* P;
		/** Matrix holding linearities */
		lrs_mp_matrix Lin;
	}; /* class lrs */
	
} /* namespace lrs */

#endif /* _LRS_HPP_ */
