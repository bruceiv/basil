#ifndef _GRAM_HPP_
#define _GRAM_HPP_

#include <functional>
#include <ostream>

#include <boost/shared_ptr.hpp>

#include "lrs/cobasis.hpp"
#include "lrs/matrix.hpp"

namespace basil {
	
	/** Matrix for gram calculations. This will be an nxn square matrix of 
	 *  integer representations of the angles between vectors of the input 
	 *  matrix. The integers chosen to represent the angles are guaranteed to 
	 *  be in one-to-one correspondance with the unique angles, but no further 
	 *  guarantees are made.
	 */
	class gram_matrix {
	friend class gram_matrix_hash;
	public:
		/** Constructs a zeroed gram matrix of dimension nxn. 
		 *  @param n		The dimension of the gram matrix
		 */
		gram_matrix(long n = 0);
		
		/** Copy constructor */
		gram_matrix(gram_matrix const& that);
		
		/** Destructor */
		~gram_matrix();
		
		/** Assignment operator */
		gram_matrix& operator= (gram_matrix const& that);
		
		/** Indexing operator.
		 *  @param i		The row index to retrieve
		 *  @param j		the column index to retrieve
		 *  @return a reference to the element at (i,j)
		 */
		int& operator() (long i, long j);
		
		/** Indexing operator.
		 *  @param i		The row index to retrieve
		 *  @param j		the column index to retrieve
		 *  @return The element at (i,j)
		 */
		int operator() (long i, long j) const;
		
		/** Takes the absolute value of every value in this matrix. This is 
		 *  equivalent to asserting that for each m_i in the generating matrix, 
		 *  m_i is equivalent to -m_i (true for arrangements, false for 
		 *  polytopes).
		 *  @return a reference to this matrix (now non-negative)
		 */
		gram_matrix& abs();
		
		/** @return the dimension of the gram matrix */
		long dim() const;
		
		/** Computes the restriction of the matrix to a given set of row and 
		 *  column indices.
		 *  @param s		The set of indices to restrict the matrix to. The 
		 *  				maximum index in s should be less than or equal to 
		 * 					n
		 *  @return a matrix R such that R[i][j] = this[s[i],s[j]]
		 */
		gram_matrix restriction(lrs::index_set s) const;
		
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
		/** The dimension of the gram matrix */
		long n;
		/** The matrix */
		int **m;
		/** Backing storage for the matrix (to save allocations) */
		int *m_;
	}; /* class gram_matrix */
	
	typedef boost::shared_ptr<gram_matrix> gram_matrix_ptr;
	
	/** Constructs the gram matrix for a given input matrix.
	 *  @param m			The input matrix, which must have no rows which are 
	 *  					the zero vector. If a zero vector is supplied, 
	 *  					undefined results ensue.
	 *  @param factorize	Perform prime factorization for an exact answer. 
	 *  					This can be set to false for quicker gram matrix 
	 *  					generation, but may give incorrect answers [true].
	 *  @return a gram matrix for the given input matrix
	 */
	gram_matrix constructGram(lrs::matrix_mpq const& m, bool factorize = true);
	
	/** Functional to hash a gram_matrix. */
	class gram_matrix_hash 
			: public std::unary_function<gram_matrix, std::size_t> {
	public:
		std::size_t operator() (gram_matrix const& m) const;
	}; /* class gram_matrix_hash */
	
} /* namespace basil */

#endif /* _GRAM_HPP_ */