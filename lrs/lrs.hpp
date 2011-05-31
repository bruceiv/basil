#ifndef _LRS_HPP_
#define _LRS_HPP_

#include "lrslib.h"
#define LRSXX_MP_H "lrsgmp.h"
#include LRSXX_MP_H
// because the LRS author saw fit to polute the global namespace with this ....
#undef copy

#include <istream>
#include <ostream>

#include <gmpxx.h>

/** Namespace for C++ LRS wrapper */
namespace lrs {
	
	/* Bring these types into the wrapper namespace to C++-ize them */
	
	/** LRS lib scalar type. */
	typedef
		lrs_mp_t
		val_t;
	/** LRS lib vector type.
	 *  Application of the indexing operator yeilds a val_t.
	 */
	typedef 
		lrs_mp_vector
		vector_t;
	/** LRS lib matrix type.
	 *  Application of the indexing operator yeilds a vector_t.
	 */
	typedef
		lrs_mp_matrix
		matrix_t;
	/** Index type for LRS vector / matrix */
	typedef
		long
		ind;
	
	/** Differentiates between expressions of equality, and expressions of 
	 *  equations */
	enum exp_type { eq = 0L, ge = 1L };
	
	/** Wraps an LRS-compatible matrix.
	 */
	class matrix {
	friend class matrix_row;
	friend class matrix_ind;
	
	public:
		/** Proxy class for a matrix index (two-dimensional) */
		class matrix_ind {
			friend class matrix;
			friend class matrix_row;
			friend class const_matrix_ind;
			
		public:
			/** Cast to GMP rational.
			 *  Copies the value in the underlying matrix.
			 */
			operator mpq_class () const;
			
			/** Assignment from GMP rational.
			 *  This modifies the underlying matrix.
			 */
			matrix_ind& operator= (mpq_class const&);
			
			/** input operator */
			friend std::istream& operator>> (std::istream& in, 
											 matrix_ind& m_i);
			
			/** output operator */
			friend std::ostream& operator<< (std::ostream& out, 
											 matrix_ind const& m_i);
			
		private:
			matrix_ind(matrix* m, ind i, ind j) : m(m), i(i), j(j) {}
			
			/** Matrix this row belongs to */
			matrix* const m;
			/** Row index in the matrix */
			ind const i;
			/** Column index in the matrix */
			ind const j;
		};
		
		/** Proxy class for a matrix index (two-dimensional) */
		class const_matrix_ind {
			friend class matrix;
			friend class matrix_row;
			friend class const_matrix_row;
			
		public:
			/** Cast to GMP rational.
			 *  Copies the value in the underlying matrix.
			 */
			operator mpq_class () const;
			
			/** output operator */
			friend std::ostream& operator<< (std::ostream& out, 
											 const_matrix_ind const& m_i);
			
		private:
			const_matrix_ind(matrix const * const m, ind i, ind j) 
				: m(m), i(i), j(j) {}
			const_matrix_ind(matrix_ind const & m_i) 
				: m(m_i.m), i(m_i.i), j(m_i.j) {}
			
			/** Matrix this row belongs to */
			matrix const * const m;
			/** Row index in the matrix */
			ind const i;
			/** Column index in the matrix */
			ind const j;
		};
		
		/** Proxy class for matrix row */
		class matrix_row {
			friend class matrix;
			friend class const_matrix_row;
			
		public:
			/** Indexing operator */
			matrix_ind operator[] (ind j) { 
				return matrix_ind(m, i, j);
			}
			const_matrix_ind const operator[] (ind j) const {
				return const_matrix_ind(m, i, j);
			}
			
			/** @return the numerator vector of the underlying matrix */
			vector_t& num();
			/** @return the denominater vector of the underlying matrix */
			vector_t& den();
			
		private:
			matrix_row(matrix* const m, ind i) : m(m), i(i) {}
			
			/** Matrix this row belongs to */
			matrix* const m;
			/** Row index in the matrix */
			ind const i;
		};
		
		/** Proxy class for constant matrix row */
		class const_matrix_row {
			friend class matrix;
			
		public:
			/** Indexing operator */
			const_matrix_ind const operator[] (ind j) const {
				return const_matrix_ind(m, i, j);
			}
			
		private:
			const_matrix_row(const matrix* const m, ind i) : m(m), i(i) {}
			
			const_matrix_row(matrix_row const & m_r) : m(m_r.m), i(m_r.i) {}
			
			/** Matrix this row belongs to */
			const matrix* const m;
			/** Row index in the matrix */
			ind const i;
		};
		
		/** Constructor.
		 *  @param n	the number of rows in the matrix
		 *  @param d	the number of columns in the matrix
		 */
		matrix(ind n, ind d);
		
		/** Destructor */
		~matrix();
		
		/* TODO put in copy constructor, assignment operator */
		
		/** Indexing operator */
		matrix_row operator[] (ind i) { 
			return matrix_row(this, i);
		}
		const_matrix_row const operator[] (ind i) const {
			return const_matrix_row(this, i);
		}
		
		/** output operator */
		friend std::ostream& operator<< (std::ostream& out, matrix const & m);
		
		/** number of rows in matrix */
		ind const n;
		/** number of columns in matrix */
		ind const d;
		
	private:
		
		/** the matrix data: numerators */
		matrix_t const& num;
		/** the matrix data: denominators */
		matrix_t const& den;
	};
	
	/** C++ wrapper class for the LRS library. */
	class lrs {
	public:
		/** Constructor / initializer.
		 *  Returns null if cannot be properly initialized.
		 */
		lrs();
		
		/** destructor */
		~lrs();
		
		/** Loads a matrix into LRS.
		 *  @param m		The matrix to load (passed by reference)
		 *  @returns true for success, false for failure
		 */
		bool loadMatrix(const matrix& m);
		
	private:
		/** Structure for holding static problem data */
		lrs_dat* Q;
		/** Structure for holding current dictionary and indices */
		lrs_dic* P;
	}; /* class lrs */
	
} /* namespace lrs */

#endif /* _LRS_HPP_ */
