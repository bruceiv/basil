#ifndef _BASIL_METRIC_HPP_
#define _BASIL_METRIC_HPP_

/** Distance metrics for d-dimensional real space, useful for computing Gram
 *  matrices.
 *
 *  @author Aaron Moss
 */

#include "basil.hpp"
#include "lrs/matrix.hpp"

namespace basil {

	/** Computes the inner product matrix of a matrix.
	 *  @param M		The matrix to compute the inner product of
	 *  @return a matrix P such that P[i][j] = inner_prod(this[i], this[j])
	 */
	lrs::matrix_mpq inner_prod_mat(lrs::matrix_mpq const& M);

	/** Computes the Q-matrix for a matrix.
	 *  @param M		The matrix to compute the Q-matrix of
	 *  @return a matrix Q such that
	 *  		Q = sum 1 to n of outer_prod(this[i], this[j])
	 */
	lrs::matrix_mpq q_mat(lrs::matrix_mpq const& M);

	/** Computes the orthogonal augmentation of a matrix. The augmentation
	 *  matrix will have full rank, and all the added rows will be orthogonal
	 *  to the row vectors already in the matrix.
	 *  @param M			The matrix to augment
	 *  @param augSigned	Should the added augment vectors be considered to
	 *   					be signed. If so, they will be paired with their
	 *   					negations (should be true for polyhedra, false for
	 *   					hyperplane arrangements) [default false]
	 *  @return a matrix augmented as above
	 */
	lrs::matrix_mpq ortho_augment(lrs::matrix_mpq const& M,
								  bool augSigned = false);

} /* namespace basil */

#endif /* BASIL_METRIC_HPP */
