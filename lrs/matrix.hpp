#ifndef _LRS_MATRIX_HPP_
#define _LRS_MATRIX_HPP_

#include <cstdlib>
#include <istream>
#include <new>
#include <ostream>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <gmpxx.h>

#include "clrs.hpp"
#include "cobasis.hpp"

namespace lrs {
	
	/** Wraps an LRS-compatible multi-precision integer vector. */
	class vector_mpz {
	//assign these friends so they can get at the v member
	friend class lrs;
	public:
		/* typedef's for STL compatibility */
		typedef val_t&			reference;
		typedef val_t const&	const_reference;
		typedef val_t* 			iterator;
		typedef val_t const*	const_iterator;
		typedef ind				size_type;
		typedef ind				difference_type;
		typedef val_t			value_type;
		typedef val_t*			pointer;
		typedef val_t const*	const_pointer;
		
		/** Constructs a vector with the given dimension. 
		 *  @param d		The dimension of the vector
		 */
		vector_mpz(size_type d);
		
		/** Copy constructor.
		 *  @param that		The vector to copy
		 */
		vector_mpz(vector_mpz const& that);
		
		/** Destructs a vector with the given dimension. */
		~vector_mpz();
		
		/** Assignment operator.
		 *  @param that		The vector to copy into this one.
		 */
		vector_mpz& operator= (vector_mpz const& that);
		
		/** Gets an iterator at the first element of the vector. */
		iterator begin();
		
		/** Gets a constant iterator at the first element of the vector. */
		const_iterator begin() const;
		
		/** Gets an iterator one past the last element of the vector. */
		iterator end();
		
		/** Gets a constant iterator one past the last element of the vector. */
		const_iterator end() const;
		
		/** Indexing operator. Returns a mutable element reference.
		 *  @param i		The index of the element to return.
		 */
		reference operator[] (size_type i);
		
		/** Indexing operator. Returns a constant element reference.
		 *  @param i		The index of the element to return.
		 */
		const_reference operator[] (size_type i) const;
		
		/** Returns the vector obtained from dividing every element in v by s.
		 *  Uses truncating integer arithmetic.
		 *  @param v		The vector to divide
		 *  @param s		The value to divide it by
		 */
		friend vector_mpz operator/ (vector_mpz const& v, const_reference s);
		
		/** Prints the vector with space-separated elements inside square 
		 *  brackets.
		 *  @param o		The output stream to print on
		 *  @param v		The vector to print
		 *  @return the same output stream
		 */
		friend std::ostream& operator<< (std::ostream& o, vector_mpz const& v);
		
		/** Standard comparison operators. Uses lexicographic comparison to 
		 *  compare the first vector to the second.
		 */
		friend bool operator<  (vector_mpz const& a, vector_mpz const& b);
		friend bool operator== (vector_mpz const& a, vector_mpz const& b);
		friend bool operator>  (vector_mpz const& a, vector_mpz const& b);
		friend bool operator<= (vector_mpz const& a, vector_mpz const& b);
		friend bool operator!= (vector_mpz const& a, vector_mpz const& b);
		friend bool operator>= (vector_mpz const& a, vector_mpz const& b);
		
		/** Returns the normalization of the given vector. The normalization is 
		 *  defined to be the vector divided by its first non-zero coordinate.
		 */
		vector_mpz normalization() const;
		
	private:
		
		/** Compares this vector to another. Uses a lexicographic comparison.
		 *  @param that		The vector to compare to
		 *  @return -1 for this vector less than that, eq 0 for this vector 
		 * 		equal to that, or gt 1 for this vector greater than that.
		 */
		int compare(vector_mpz const& that) const;
		
		/** Internal storage of vector data */
		vector_t v;
		/** Dimension of the vector. */
		size_type d;
	};
	
	/** Functional to hash a vector_mpz */
	class vector_mpz_hash 
			: public std::unary_function<vector_mpz, std::size_t> {
	public:
		
		std::size_t operator() (vector_mpz const& v) const {
			std::size_t seed = 0UL;
			vector_mpz& v_n = const_cast<vector_mpz&>(v);
			for (vector_mpz::iterator i = v_n.begin(); i != v_n.end(); ++i) {
				boost::hash_combine(seed, mptoi(*i) );
			}
			return seed;
			/*
			return boost::hash_range(
				boost::make_transform_iterator( 
					v.begin(), 
					boost::pointer_to_unary_function<val_t const, long>(
						getLong) 
				),
				boost::make_transform_iterator( 
					v.end(), 
					boost::pointer_to_unary_function<val_t const, long>(
						getLong) 
				)
			);
			*/
		}
	
	};
	
	/** Wraps an LRS-compatible matrix. */
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
		 *  @param n			the number of rows in the matrix
		 *  @param d			the number of columns in the matrix
		 *  @throw bad_alloc	if allocation of memory for the matrix fails
		 */
		matrix(ind n, ind d) throw(std::bad_alloc);
		
		/** Copy constructor. 
		 *  @param that			matrix to copy
		 *  @throw bad_alloc	if allocation of memory for the matrix fails
		 */
		matrix(matrix const& that) throw(std::bad_alloc);
		
		/** Destructor */
		~matrix();
		
		/** Assignment operator.
		 *  @param that			matrix to copy into this one
		 *  @throw bad_alloc	if allocation of memory for the matrix fails
		 */
		matrix& operator= (matrix const& that) throw(std::bad_alloc);
		
		/** Indexing operator */
		matrix_row operator[] (ind i) { 
			return matrix_row(this, i);
		}
		const_matrix_row const operator[] (ind i) const {
			return const_matrix_row(this, i);
		}
		
		/** output operator */
		friend std::ostream& operator<< (std::ostream& out, matrix const & m);
		
		/** Size accessor: rows */
		ind const & n() const { return n_; }
		
		/** Size accessor: columns */
		ind const & d() const { return d_; }
		
	private:
		/** number of rows in matrix */
		ind n_;
		/** number of columns in matrix */
		ind d_;
		
		/** the matrix data: numerators */
		matrix_t const * num;
		/** the matrix data: denominators */
		matrix_t const * den;
	};
	
} /* namespace lrs */

#endif /* _LRS_MATRIX_HPP_ */