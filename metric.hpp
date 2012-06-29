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
	 *  common factors, and d and r are both positive.
	 */
	struct mpr {

		/** Constructor; initializes the mpr to the canonical form of 0
		 *  [0*sqrt(1)/1] */
		mpr();

		/** Constructor; initializes mpr to n_*sqrt(r_)/d_. Performs no
		 *  normalization.
		 *  @param n		The numerator
		 *  @param r		The radical (should be positive)
		 *  @param d		The denominator (should be non-zero)
		 */
		mpr(mpz_class n_, mpz_class r_, mpz_class d_);

		/** Creates a normalized mpr.
		 *  @param n_		The numerator
		 *  @param rf		The factors of the radical
		 *  @param d_		The denominator (should be non-zero)
		 *  @param factor	The prime factorizer to use to factor r [defaults
		 *  				to a new instance - it is strongly recommended that
		 *  				the prime factorization instance be saved if more
		 *  				than one normalization is desired].
		 *  @return an mpr reprsenting the parameters, but with no factors of r
		 *  		which are perfect squares, and n/d in lowest terms with a
		 *  		positive d.
		 */
		static mpr makeNorm(mpz_class n_, prime::factor_list rf, mpz_class d_,
						  prime::factorizer factor = prime::factorizer());

		/** Constructor; initializes mpr to x*sqrt(1)/1. */
		mpr(int x);

		/** Assignment operator; sets this to x*sqrt(1)/1 */
		mpr& operator= (int x);

		/** Normalizes the fraction n/d to lowest terms, with d positive. Does
		 *  not normalize the radical part r.
		 */
		void normRational();

		/** Normalizes the mpr. After this is called, r will have no factors
		 *  which are perfect squares, and n/d will be in lowest terms with a
		 *  positive d.
		 *  @param factor	The prime factorizer to use to factor r [defaults
		 *  				to a new instance - it is strongly recommended that
		 *  				the prime factorization instance be saved if more
		 *  				than one normalization is desired].
		 */
		void norm(prime::factorizer factor = prime::factorizer());

		/** Numerator */
		mpz_class n;
		/** Radical */
		mpz_class r;
		/** Denominator */
		mpz_class d;
	}; /* struct mpr */

	/** Equality operator for mpr type
	 *  @param a				the first mpr to check (should be normalized)
	 *  @param b				the second mpr to check (should be normalized)
	 *  @return true for a == b, false otherwise
	 */
	bool operator== (mpr const& a, mpr const& b);

	/** Inequality operator for mpr type
	 *  @param a				the first mpr to check (should be normalized)
	 *  @param b				the second mpr to check (should be normalized)
	 *  @return true for a != b, false otherwise
	 */
	bool operator!= (mpr const& a, mpr const& b);

	/** Output operator for mpr type */
	std::ostream& operator<< (std::ostream& o, mpr x);

	/** Takes the absolute value of an mpr
	 *  @param x				the mpr to take the absolute value of - should
	 *  						be in normal form
	 *  @return absolute value of x
	 */
	mpr abs(mpr const& x);

	/** Takes the sign of an mpr
	 *  @param x				the mpr to take the sign of - should be in
	 *  						normal form
	 *  @return -1 for x < 0, 0 for x == 0, 1 for x > 0
	 */
	int sgn(mpr const& x);

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


	/** Adds a [1 0 0 ... 0] row to a matrix. When interpreted as an LRS-style
	 *  constraint, equivalent to 1 >= 0 (always true), when interpreted as a
	 *  vector in d-space, ensures that any automorphisms will fix the x_0 = 1
	 *  plane.
	 *  @param M		The matrix to append to
	 *  @return a copy of the matrix with an extra row [1 0 0 ... 0] added
	 */
	lrs::matrix_mpq fixPlane(lrs::matrix_mpq const& M);

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
	 *   					hyperplane arrangements) [default true]
	 *  @return a matrix augmented as above
	 */
	lrs::matrix_mpq orthoAugment(lrs::matrix_mpq const& M,
								 bool augSigned = true);

	/** Computes the column rank augmentation of a matrix.
	 *  @param M			The matrix to augment
	 *  @param rows			The rows to keep and augment
	 *  @return a matrix augmented as above
	 */
	lrs::matrix_mpq colRankAugment(lrs::matrix_mpq const& M,
								   lrs::index_set const& rows);

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
	 *  @param M			The matrix to compute the inner products for.
	 *  					Should not have any zero rows - inclusion of a zero
	 *  					row will result in undefined results.
	 *  @return a matrix P such that P[i][j] =
	 *  		inner_prod(M[i],M[j])/( ||M[i]|| * ||M[j]|| )
	 */
	matrix_mpr normedInnerProdMat(lrs::matrix_mpq const& M);

	/** Selects the submatrix of a matrix determined by the given list of
	 *  indices.
	 *  @param M		The matrix to select from
	 *  @param l		The set of row indices to select from the matrix.
	 *  				The maximum index in l should be less than or equal
	 *  				to M.n.
	 *  @return a matrix R such that R[i][j] = M[l[i],j]
	 */
	lrs::matrix_mpq select_rows(lrs::matrix_mpq const& M, index_list const& l);

} /* namespace basil */

#endif /* BASIL_METRIC_HPP */
