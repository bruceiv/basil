#ifndef _FUND_DOMAIN_HPP_
#define _FUND_DOMAIN_HPP_

#include <vector>

#include "basil.hpp"
#include "lrs/matrix.hpp"

namespace basil {

	/** Represents a fundamental domain, a polyhedon which tiles a polyhedron
	 *  or arrangement under the action of some symmetry group.
	 */
	class fund_domain {
	public:
		/** Vector type */
		typedef lrs::vector_mpq vector_mpq;
		/** Matrix type */
		typedef lrs::matrix_mpq matrix_mpq;

		/** Default constructor.
		 *  @param qInv		The Q-matrix inverse (used to compute norms)
		 */
		fund_domain(matrix_mpq qInv);

		/** Adds a new constraint to the fundamental domain. Adds a halfspace h
		 *  whose bounding hyperplane is the normal bisector of the line
		 *  segment ab, such that a is included in the halfspace, but b is not.
		 *  @param a		The point to include (should be length d+1, with
		 *  				the leading coefficient 1 and all others defining
		 *  				the point in d-space in the normal order)
		 *  @param b		The point to exclude (same constraints as a)
		 */
		void add_constraint(vector_mpq const& a, vector_mpq const& b);

		/** Checks if this fundamental domain contains a given point.
		 *  @param x		The point to check, should have length d
		 *  @return true for x contained within the polyhedron, false otherwise
		 */
		bool contains(vector_mpq const& x) const;

		/** Gets the dimension of the fundamental domain's underlying space */
		ind dim() const;
		/** Gets the size (number of halfspaces) of the fundamental domain */
		ind size() const;

	private:
		/** The polytope underlying the fundamental domain */
		std::vector<vector_mpq> p;
		/** The Q-matrix inverse (used to compute norms) */
		matrix_mpq qInv;
	}; /* class fund_domain */

	/** Gets a halfspace to include a but not b in a fundamental domain.
	 *  @param a		The point to include
	 *  @param b		The point to exclude
	 *  @return the halfspace including a but not b, such that the bounding
	 *  		hyperplane is the normal bisector of the line segment between
	 *  		a and b
	 */
	fund_domain::vector_mpq splitting_halfspace(
			fund_domain::vector_mpq const& a, fund_domain::vector_mpq const& b);

} /* namespace basil */

#endif /* _FUND_DOMAIN_HPP_ */
