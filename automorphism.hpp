#ifndef _AUTOMORPHISM_HPP_
#define _AUTOMORPHISM_HPP_

/** Matrix automorphism calculations, for use in determining the restricted
 *  automorphism groups of polytopes and hyperplane arrangements.
 *
 *  @author Aaron Moss
 */

/*  Copyright: Aaron Moss, 2012, moss.aaron@unb.ca  */

/*  This file is part of Basil.

    Basil is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Basil is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Basil.  If not, see <http://www.gnu.org/licenses/>.  */

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
