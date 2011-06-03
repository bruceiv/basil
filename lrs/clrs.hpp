#ifndef _CLRS_HPP_
#define _CLRS_HPP_

#include "lrslib.h"
#include "lrsgmp.h"
// because the LRS author saw fit to polute the global namespace with this ....
#undef copy

//need to declare these again, because they don't link properly otherwise
lrs_mp_matrix lrs_alloc_mp_matrix(long,long);
void lrs_clear_mp_matrix(lrs_mp_matrix,long,long);


/* Bring these types into the wrapper namespace to C++-ize them */
namespace lrs {
	
	/** Copies src to dst. 
	 *  Replaces the copy macro in lrsgmp.h that I undef'd from the global 
	 *  namespace.
	 */
	static void copy(mpz_t dst, mpz_t src) { mpz_set(dst, src); }
	
	/** LRS lib scalar type. */
	typedef
		lrs_mp_t
		val_t;
	/** LRS lib vector type.
	 *  Application of the indexing operator yeilds a val_t.
	 */
	typedef 
		lrs_mp_vector
		vector_t;
	/** LRS lib matrix type.
	 *  Application of the indexing operator yeilds a vector_t.
	 */
	typedef
		lrs_mp_matrix
		matrix_t;
	/** Index type for LRS vector / matrix. */
	typedef
		long
		ind;
}

#endif /* _CLRS_HPP_ */