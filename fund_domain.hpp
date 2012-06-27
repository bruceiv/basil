#ifndef _FUND_DOMAIN_HPP_
#define _FUND_DOMAIN_HPP_

/** Calculations on fundamental domains.
 *
 *  @author Aaron Moss
 */

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
		/** Iterator type */
		typedef std::vector<vector_mpq>::const_iterator const_iterator;

		/** Empty constructor.
		 *  Creates a fundamental domain with an empty Q-matrixs
		 */
		fund_domain();

		/** Default constructor.
		 *  @param qInv		The Q-matrix inverse (used to compute norms)
		 */
		fund_domain(matrix_mpq qInv);

		/** Builds the fundamental domain from a seed vertex and the group.
		 *  Adds the constraints generated from the seed vertex and its images
		 *  under some list of permutations to the fundamental domain.
		 *  @param s		The seed vertex
		 *  @param sBasis	The set of basis rows for the seed
		 *  @param A		The constraint matrix the seed is in the context of
		 *  @param l		The list of permutations to apply to the seed
		 */
		void build_from_seed(vector_mpq const& s, index_set const& sBasis,
				matrix_mpq const& A, permutation_list const& l);

		/** Creates a new constraint to the fundamental domain.
		 *  @param a		The point to include (should be length d+1, with
		 *  				the leading coefficient 1 and all others defining
		 *  				the point in d-space in the normal order)
		 *  @param b		The point to exclude (same constraints as a)
		 *  @param c		A halfspace whose bounding hyperplane is the normal
		 *  				bisector of the line segment ab, such that a is
		 *  				included in the halfspace, but b is not.
		 */
		vector_mpq get_constraint(vector_mpq const& a,
				vector_mpq const& b) const;

		/** Adds a new constraint to the fundamental domain.
		 *  @param c		The constraint to add
		 */
		void push_back(vector_mpq const& c);

		/** Checks if this fundamental domain contains a given point.
		 *  @param x		The point to check, should have length d
		 *  @return true for x contained within the polyhedron, false otherwise
		 */
		bool contains(vector_mpq const& x) const;

		/** @return the constraints stored in this fundamental domain */
		std::vector<vector_mpq> const& constraints() const;

		/** @return an iterator pointing to the first element in the
		 *  fundamental domain */
		const_iterator begin() const;

		/** @return an iterator pointing beyond the last element in the
		 *  fundamental domain */
		const_iterator end() const;

		/** Gets the dimension of the fundamental domain's underlying space */
		uind dim() const;
		/** Gets the size (number of halfspaces) of the fundamental domain */
		uind size() const;

	private:
		/** The polytope underlying the fundamental domain */
		std::vector<vector_mpq> p;
		/** The Q-matrix inverse (used to compute norms) */
		matrix_mpq qInv;
	}; /* class fund_domain */

} /* namespace basil */

#endif /* _FUND_DOMAIN_HPP_ */
