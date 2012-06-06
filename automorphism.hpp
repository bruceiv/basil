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
	
	/** Takes a group of degree m and makes a group of degree n from it.
	 *  @param g			The group to shrink
	 *  @param n			The degree of the shrunken group
	 *  @return a group expressed as a subgroup of S_n, excluding any elements
	 *  		of the original group that do not setwise fix S_n.
	 */
	permutation_group_ptr shrink_group_to(permutation_group const& g, uind n);

} /* namespace basil */

#endif /* _AUTOMORPHISM_HPP_ */
