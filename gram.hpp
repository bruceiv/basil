#ifndef _GRAM_HPP_
#define _GRAM_HPP_

/** Computation and manipulation of "Gram" (normed inner product
 *  representative) matrices for the purpose of symmetry calculations.
 *
 *  @author Aaron Moss
 */

#include <functional>
#include <ostream>

#include <boost/shared_ptr.hpp>

#include "basil.hpp"
#include "metric.hpp"

namespace basil {
	
	/** Matrix for Gram calculations. This will be an nxn square matrix of
	 *  non-negative integer representations of the angles between vectors of 
	 *  the input matrix. The integers chosen to represent the angles are 
	 *  guaranteed to be in one-to-one correspondance with the unique angles, 
	 *  and all values i, 0 <= i < k, are guaranteed to be included in any 
	 *  unrestricted Gram matrix, but no further guarantees are made.
	 */
	class gram_matrix {
	friend class gram_matrix_hash;
	public:
		/** Constructs a zeroed gram matrix of dimension nxn. 
		 *  @param n		The dimension of the Gram matrix
		 *  @param k_		The maximum value in the Gram matrix
		 */
		gram_matrix(uind n = 0, uind k_ = 0);
		
		/** Copy constructor */
		gram_matrix(gram_matrix const& that);
		
		/** Destructor */
		~gram_matrix();
		
		/** Assignment operator */
		gram_matrix& operator= (gram_matrix const& that);
		
		/** @return k, the maximum value in this matrix. Included for 
		 *  compatibility with PermLib */
		uind& k() { return k_; }
		
		/** @return k, the maximum value in this matrix. Included for 
		 *  compatibility with PermLib */
		uind k() const { return k_; }
		
		/** Indexing operator.
		 *  @param i		The row index to retrieve
		 *  @param j		the column index to retrieve
		 *  @return a reference to the element at (i,j)
		 */
		inline int& operator() (uind i, uind j)
			{ return at(i, j); }
		
		/** Indexing operator.
		 *  @param i		The row index to retrieve
		 *  @param j		the column index to retrieve
		 *  @return The element at (i,j)
		 */
		int operator() (uind i, uind j) const
			{ return at(i, j); }
		
		/** Synonym for operator(), for compatibility with PermLib */
		int& at(unsigned long i, unsigned long j);
		
		/** Synonym for operator(), for compatibility with PermLib */
		int at(unsigned long i, unsigned long j) const;
		
		/** @return the dimension of the gram matrix */
		uind dim() const;
		
		/** Synonym for dim(), for compatibility with PermLib */
		inline unsigned long dimension() const { return dim(); }
		
		/** Computes the restriction of the matrix to a given set of row and 
		 *  column indices.
		 *  @param s		The set of indices to restrict the matrix to. The 
		 *  				maximum index in s should be less than or equal to 
		 * 					n
		 *  @return a matrix R such that R[i][j] = this[s[i],s[j]]
		 */
		gram_matrix restriction(index_set s) const;
		
		/** @return the Gram matrix equivalent to taking the absolute value of
		 *  every value in this matrix.
		 */
		gram_matrix abs() const;
		
		/** @return the Gram matrix where each row is doubled with it's
		 *  negation. (That is, will return a Gram matrix D = (d_i,j) such for
		 *  the original matrix A = (a_i,j), a_i,j = x => d_2i,2j = x, 
		 *  d_2i+1,2j+1 = x, d_2i,2j+1 = -x, d_2i+1,2j = -x ) 
		 */
		gram_matrix doubled() const;
		
		/** @return the Gram matrix in a format that PermLib can deal with it
		 *  in, namely with all values in the matrix in the range [0,k), and 
		 *  each value in that range appearing at least once.
		 */
		gram_matrix permlibCanon() const;
		
		/** Sorts this matrix, first sorting each row, then lexicographically 
		 *  sorting the matrix by rows. This has the effect of providing a 
		 *  canonical permuation of a matrix, such that two matrices restricted 
		 *  from the same base matrix can be compared for equality.
		 *  @return a reference to this matrix (now sorted)
		 */
		gram_matrix& sort();
		
		/** Tests two gram matrices for equality, elementwise */
		friend bool operator== (gram_matrix const& a, gram_matrix const& b);
		
		/** Prints a gram matrix */
		friend std::ostream& operator<< (std::ostream& o, gram_matrix const& g);
		
	private:
		/** The dimension of the Gram matrix */
		uind n;
		/** The maximum value in the range of unique values in this matrix */
		uind k_;
		/** The matrix */
		int **m;
		/** Backing storage for the matrix (to save allocations) */
		int *m_;
	}; /* class gram_matrix */
	
	typedef boost::shared_ptr<gram_matrix> gram_matrix_ptr;
	
	/** Constructs a Gram matrix from the given rational matrix.
	 *  @param M			The matrix of angles to represent with integer
	 *  					representations - should be square.
	 *  @return A Gram matrix with a unique integer representation for each
	 *  		value. If M[i,j] is represented by x, -M[i,j] will be
	 *  		represented by -x.
	 */
	gram_matrix constructGram(lrs::matrix_mpq const& M);

	/** Constructs a Gram matrix from the given radical matrix.
	 *  @param M			The matrix of angles to represent with integer
	 *  					representations - should be square.
	 *  @return A Gram matrix with a unique integer representation for each
	 *  		value. If M[i,j] is represented by x, -M[i,j] will be
	 *  		represented by -x.
	 */
	gram_matrix constructGram(matrix_mpr const& M);

	/** Functional to hash a gram_matrix. */
	class gram_matrix_hash 
			: public std::unary_function<gram_matrix, std::size_t> {
	public:
		std::size_t operator() (gram_matrix const& m) const;
	}; /* class gram_matrix_hash */
	
} /* namespace basil */

#endif /* _GRAM_HPP_ */
