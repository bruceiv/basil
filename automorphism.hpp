#ifndef _AUTOMORPHISM_HPP_
#define _AUTOMORPHISM_HPP_

/** Matrix automorphism calculations, for use in determining the restricted
 *  automorphism groups of polytopes and hyperplane arrangements.
 *
 *  @author Aaron Moss
 */

#include "basil.hpp"
#include "gram.hpp"

namespace basil {
	
	/** Computes the restricted automorphism group for the polytope with the 
	 *  given gram matrix.
	 *  @param g			The gram matrix of the polytope to compute the 
	 *  					restricted automorphism group for.
	 *  @return the restricted automorphism group for the polytope
	 */
	permutation_group_ptr compute_restricted_automorphisms(
		gram_matrix const& g);
	
	/** Computes the restricted automorphism group for the arrangement with the 
	 *  given gram matrix.
	 *  @param g			The gram matrix of the arrangement to compute the 
	 *  					restricted automorphism group for.
	 *  @return the restricted automorphism group for the arrangement
	 */
	permutation_group_ptr compute_arrangement_automorphisms(
		gram_matrix const& g);
	
} /* namespace basil */

#endif /* _AUTOMORPHISM_HPP_ */
