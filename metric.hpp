#ifndef _BASIL_METRIC_HPP_
#define _BASIL_METRIC_HPP_

/** Distance metrics for d-dimensional real space, useful for computing Gram
 *  matrices.
 *
 *  @author Aaron Moss
 */

#include <gmpxx.h>

#include "basil.hpp"
#include "lrs/matrix.hpp"
#include "prime/factors.hpp"

namespace basil {

	/** Represents a multiprecision radical fraction (a number of the form
	 *  n*sqrt(r)/d). In normalized form, the radical r has no factors which
	 *  are perfect squares, and the numerator and denominator n and d have no
	 *  common factors.
	 */
	struct mpr {

		/** Constructor; initializes the mpr to the canonical form of 0
		 *  [0*sqrt(1)/1] */
		mpr();

		/** Constructor; initializes mpr to n_*sqrt(r_)/d_. Performs no
		 *  normalization. */
		mpr(mpz_class n_, mpz_class r_, mpz_class d_);

		/** Normalizes the fraction n/d to lowest terms. Does not normalize the
		 *  radical part r.
		 */
		void normRational();



		/** Numerator */
		mpz_class n;
		/** Radical */
		mpz_class r;
		/** Denominator */
		mpz_class d;
	}; /* struct mpr */

	/** Equality operator for mpr type */
	bool operator== (mpr const& a, mpr const& b);

	/** Inequality operator for mpr type */
	bool operator!= (mpr const& a, mpr const& b);

	/** Output operator for mpr type */
	std::ostream& operator<< (std::ostream& o, mpr x);

	/** Matrix of mpr objects. Designed for minimal compatibility with
	 *  lrs::matrix_mpq.
	 */
	class matrix_mpr {
	public:
		/** Constructs a matrix with the given number of rows and columns, with
		 *  all elements initialized to zero. If called as the default
		 *  constructor, constructs an empty matrix.
		 *  @param n			the number of rows in the matrix [default 0]
		 *  @param d			the number of columns in the matrix [default 0]
		 */
		matrix_mpr(ind n = 0, ind d = 0);

		/** Copy constructor.
		 *  @param that			matrix to copy
		 */
		matrix_mpr(matrix_mpr const& that);

		/** Destructor */
		~matrix_mpr();

		/** Assignment operator.
		 *  @param that			matrix to copy into this one
		 */
		matrix_mpr& operator= (matrix_mpr const& that);

		/** Gets the number of rows in this matrix. */
		ind size() const;

		/** Gets the dimension of the rows in this matrix. */
		ind dim() const;

		/** Gets the (i,j)th element of this matrix */
		mpr& elem(ind i, ind j);

		/** Gets the (i,j)th element of this matrix */
		mpr const& elem(ind i, ind j) const;

	private:
		/** Matrix data storage */
		mpr* m;
		/** Number of rows in the matrix */
		ind n;
		/** Dimension of the matrix rows */
		ind d;
	}; /* class matrix_mpr */

	/** Computes the inner product matrix of a matrix.
	 *  @param M		The matrix to compute the inner product of
	 *  @return a matrix P such that P[i][j] = inner_prod(M[i], M[j])
	 */
	lrs::matrix_mpq innerProdMat(lrs::matrix_mpq const& M);

	/** Computes the inverse of the Q-matrix of a matrix.
	 *  @param M		The matrix to compute the Q-matrix inverse for
	 *  @return a matrix inv(Q) such that
	 *  		Q = sum 1 to n of outer_prod(M[i], M[j])
	 */
	lrs::matrix_mpq invQMat(lrs::matrix_mpq const& M);

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
	lrs::matrix_mpq orthoAugment(lrs::matrix_mpq const& M,
								 bool augSigned = false);

	/** Computes the inner product matrix of M*T with M.
	 *  @param M			The matrix to compute the inner products for
	 *  @param T			The matrix to transform M by (should be M.d*M.d
	 *  					square)
	 *  @return a matrix P such that P[i][j] = inner_prod(M[i]*T,M[j])
	 */
	lrs::matrix_mpq transformedInnerProdMat(lrs::matrix_mpq const& M,
											lrs::matrix_mpq const& T);

	/** Computes the inner product matrix of a matrix, normalized according to
	 *  the Euclidean metric.
	 *  @param M			The matrix to compute the inner products for
	 *  @return a matrix P such that P[i][j] =
	 *  		inner_prod(M[i],M[j])/( ||M[i]|| * ||M[j]|| )
	 */
	matrix_mpr normedInnerProdMat(lrs::matrix_mpq const& M);

} /* namespace basil */

#endif /* BASIL_METRIC_HPP */
