#ifndef _LRS_MATRIX_HPP_
#define _LRS_MATRIX_HPP_

#include <cstdlib>
#include <functional>
#include <istream>
#include <new>
#include <ostream>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include "clrs.hpp"
#include "cobasis.hpp"

namespace lrs {
	
	class vector_mpz;
	class matrix_mpq;
	
	/** Wraps a multi-precision rational vector. */
	class vector_mpq {
	friend class vector_mpz;
	friend class matrix_mpq;
	public:
		/* typedef's for STL compatibility */
		typedef mpq_t&			reference;
		typedef mpq_t const&	const_reference;
		typedef mpq_t* 			iterator;
		typedef mpq_t const*	const_iterator;
		typedef ind				size_type;
		typedef ind				difference_type;
		typedef mpq_t			value_type;
		typedef mpq_t*			pointer;
		typedef mpq_t const*	const_pointer;
		
		/** Constructs a vector with the given dimension. 
		 *  @param d		The dimension of the vector
		 */
		vector_mpq(size_type d);
		
		/** Copy constructor.
		 *  @param that		The vector to copy
		 */
		vector_mpq(vector_mpq const& that);
		
		/** Copy constructor.
		 *  @param that		The integer vector to copy. Will initialize all 
		 * 					denominators to 1
		 */
		vector_mpq(vector_mpz const& that);
		
		/** Struct constructor; allows creation of a vector_mpq view of 
		 *  pre-allocated memory. A vector created this way will not deallocate 
		 *  the shared memory, though it may modify the values contained 
		 *  therein */
		vector_mpq(mpq_t* v, size_type d);
		
		/** Destructs a vector. */
		~vector_mpq();
		
		/** Assignment operator.
		 *  @param that		The vector to copy into this one. Undefined results 
		 * 					if this vector is not self-allocated and that.d 
		 * 					differs from this->d.
		 */
		vector_mpq& operator= (vector_mpq const& that);
		
		/** Assignment operator.
		 *  @param that		The integer vector to copy into this one. Will 
		 * 					initialize all denominators to 1. Undefined results 
		 * 					if this vector is not self-allocated and that.d 
		 * 					differs from this->d.
		 */
		vector_mpq& operator= (vector_mpz const& that);
		
		/** Gets an iterator at the first element of the vector. */
		iterator begin();
		
		/** Gets a constant iterator at the first element of the vector. */
		const_iterator begin() const;
		
		/** Gets an iterator one past the last element of the vector. */
		iterator end();
		
		/** Gets a constant iterator one past the last element of the vector. */
		const_iterator end() const;
		
		/** Gets the dimension of this vector. */
		size_type size() const;
		
		/** Indexing operator. Returns a mutable element reference.
		 *  @param i		The index of the element to return.
		 */
		reference operator[] (size_type i);
		
		/** Indexing operator. Returns a constant element reference.
		 *  @param i		The index of the element to return.
		 */
		const_reference operator[] (size_type i) const;
		
		/** Prints the vector with space-separated elements inside square 
		 *  brackets.
		 *  @param o		The output stream to print on
		 *  @param v		The vector to print
		 *  @return the same output stream
		 */
		friend std::ostream& operator<< (std::ostream& o, vector_mpq const& v);
		
		/** Less than. Equivalent to a.compare(b) < 0 */
		friend bool operator<  (vector_mpq const& a, vector_mpq const& b);
		/** Equal. Equivalent to a.compare(b) == 0 */
		friend bool operator== (vector_mpq const& a, vector_mpq const& b);
		/** Greater than. Equivalent to a.compare(b) > 0 */
		friend bool operator>  (vector_mpq const& a, vector_mpq const& b);
		/** Less than or equal. Equivalent to a.compare(b) <= 0 */
		friend bool operator<= (vector_mpq const& a, vector_mpq const& b);
		/** Not equal. Equivalent to a.compare(b) != 0 */
		friend bool operator!= (vector_mpq const& a, vector_mpq const& b);
		/** Greater than or equal. Equivalent to a.compare(b) >= 0 */
		friend bool operator>= (vector_mpq const& a, vector_mpq const& b);
		
		/** Computes the inner product of two vectors.
		 *  @param a		A vector of length d
		 *  @param b		A vector of length d
		 *  @return \f$\sum_{i=0}^{d} a_i * b_i\f$
		 *  @throws std::runtime_error on a.d != b.d
		 */
		friend mpq_class inner_prod (vector_mpq const& a, vector_mpq const& b);
		
		/** Gets the vector of numerators of this vector */
		vector_mpz num();
		/** Gets the vector of denominators of this vector */
		vector_mpz den();
		
	private:
		
		/** Compares this vector to another, lexicographically.
		 *  @param that		The vector to compare to
		 *  @return -1 for this vector less than that, 0 for this vector equal 
		 *  		to that, or 1 for this vector greater than that.
		 */
		int compare(vector_mpq const& that) const;
		
		/** Internal storage of vector data */
		mpq_t* v;
		/** Dimension of the vector. */
		size_type d;
		/** Is this vector responsible for its own storage? */
		bool selfAlloc;
	};
	
	/** Functional to hash a vector_mpq.  */
	class vector_mpq_hash 
			: public std::unary_function<vector_mpq, std::size_t> {
	public:
		std::size_t operator() (vector_mpq const& v) const;
	};
	
	/** Wraps an LRS-compatible multi-precision integer vector. */
	class vector_mpz {
	//assign these friends so they can get at the v member
	friend class lrs;
	friend class vector_mpq;
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
		
		/** Gets the dimension of this vector. */
		size_type size() const;
		
		/** Indexing operator. Returns a mutable element reference.
		 *  @param i		The index of the element to return.
		 */
		reference operator[] (size_type i);
		
		/** Indexing operator. Returns a constant element reference.
		 *  @param i		The index of the element to return.
		 */
		const_reference operator[] (size_type i) const;
		
		/** Prints the vector with space-separated elements inside square 
		 *  brackets.
		 *  @param o		The output stream to print on
		 *  @param v		The vector to print
		 *  @return the same output stream
		 */
		friend std::ostream& operator<< (std::ostream& o, vector_mpz const& v);
		
		/** Less than. Equivalent to a.compare(b) < 0 */
		friend bool operator<  (vector_mpz const& a, vector_mpz const& b);
		/** Equal. Equivalent to a.compare(b) == 0 */
		friend bool operator== (vector_mpz const& a, vector_mpz const& b);
		/** Greater than. Equivalent to a.compare(b) > 0 */
		friend bool operator>  (vector_mpz const& a, vector_mpz const& b);
		/** Less than or equal. Equivalent to a.compare(b) <= 0 */
		friend bool operator<= (vector_mpz const& a, vector_mpz const& b);
		/** Not equal. Equivalent to a.compare(b) != 0 */
		friend bool operator!= (vector_mpz const& a, vector_mpz const& b);
		/** Greater than or equal. Equivalent to a.compare(b) >= 0 */
		friend bool operator>= (vector_mpz const& a, vector_mpz const& b);
		
		/** Returns the rationalization of the given vector. The 
		 *  rationalization is defined to be the vector divided by its first 
		 *  coordinate (if non-zero), or the vector otherwise.
		 */
		vector_mpq rationalization() const;
		
	private:
		/** Compares this vector to another lexicographically.
		 *  @param that		The vector to compare to
		 *  @return -1 for this vector less than that, 0 for this vector equal 
		 *  		to that, or 1 for this vector greater than that.
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
		std::size_t operator() (vector_mpz const& v) const ;
	};
	
	/** Wraps a multi-precision rational matrix. */
	class matrix_mpq {
	public:
		/* typedef's for STL compatibility */
		typedef vector_mpq			reference;
		typedef vector_mpq const	const_reference;
		typedef ind					size_type;
		typedef ind					difference_type;
		typedef vector_mpq			value_type;
		typedef vector_mpq*			pointer;
		typedef vector_mpq const*	const_pointer;
		
		class const_iterator;
		/** Row iterator */
		class iterator 
				: public boost::iterator_facade<iterator, 
						matrix_mpq::value_type, 
						std::random_access_iterator_tag,
						matrix_mpq::reference,
						matrix_mpq::difference_type> {
		private:
			friend class matrix_mpq;
			friend class matrix_mpq::const_iterator;
			/* required for boost::iterator_facade */
			friend class boost::iterator_core_access;
			
			/** Constructor */
			iterator(mpq_t* m, matrix_mpq::size_type i, 
					matrix_mpq::size_type d);
			
			/** Returns the referenced vector. In compliance with 
			 *  boost::iterator_facade */
			matrix_mpq::value_type dereference() const;
			
			/** Tests the equality of two iterators. In compliance with 
			 *  boost::iterator_facade */
			bool equal(iterator& that) const;
			
			/** Tests the equality of two iterators. In compliance with 
			 *  boost::iterator_facade */
			bool equal(const_iterator& that) const;
			
			/** Moves to the next row. In compliance with 
			 *  boost::iterator_facade */
			void increment();
			
			/** Moves to the previous row. In compliance with 
			 *  boost::iterator_facade */
			void decrement();
			
			/** Moves n rows (n may be negative). In compliance with 
			 *  boost::iterator_facade */
			void advance(matrix_mpq::difference_type n);
			
			/** Returns number of rows that is ahead of this (may be negative). 
			 *  Undefined if iterators are not over the same matrix. In 
			 *  compliance with boost::iterator_facade */
			matrix_mpq::difference_type distance_to(iterator& that) const;
			
			/** Returns number of rows that is ahead of this (may be negative). 
			 *  Undefined if iterators are not over the same matrix. In 
			 *  compliance with boost::iterator_facade */
			matrix_mpq::difference_type distance_to(const_iterator& that) const;
			
			/** Pointer to the underlying matrix */
			mpq_t* m;
			/** Row index of this iterator in the matrix */
			matrix_mpq::size_type i;
			/** Dimension of the matrix rows. */
			matrix_mpq::size_type d;
		};
		/** Constant row iterator */
		class const_iterator 
				: public boost::iterator_facade<const_iterator, 
						matrix_mpq::value_type const, 
						std::random_access_iterator_tag,
						matrix_mpq::const_reference,
						matrix_mpq::difference_type> {
		private:
			friend class matrix_mpq;
			/* required for boost::iterator_facade */
			friend class boost::iterator_core_access;
			
			/** Constructor */
			const_iterator(mpq_t* m, matrix_mpq::size_type i, 
					matrix_mpq::size_type d);
			
			/** Constifying copy constructor */
			const_iterator(iterator& that);
			
			/** Returns the referenced vector. In compliance with 
			 *  boost::iterator_facade */
			matrix_mpq::value_type const dereference() const;
			
			/** Tests the equality of two iterators. In compliance with 
			 *  boost::iterator_facade */
			bool equal(iterator& that) const;
			
			/** Tests the equality of two iterators. In compliance with 
			 *  boost::iterator_facade */
			bool equal(const_iterator& that) const;
			
			/** Moves to the next row. In compliance with 
			 *  boost::iterator_facade */
			void increment();
			
			/** Moves to the previous row. In compliance with 
			 *  boost::iterator_facade */
			void decrement();
			
			/** Moves n rows (n may be negative). In compliance with 
			 *  boost::iterator_facade */
			void advance(matrix_mpq::difference_type n);
			
			/** Returns number of rows that is ahead of this (may be negative). 
			 *  Undefined if iterators are not over the same matrix. In 
			 *  compliance with boost::iterator_facade */
			matrix_mpq::difference_type distance_to(iterator& that) const;
			
			/** Returns number of rows that is ahead of this (may be negative). 
			 *  Undefined if iterators are not over the same matrix. In 
			 *  compliance with boost::iterator_facade */
			matrix_mpq::difference_type distance_to(const_iterator& that) const;
			
			/** Pointer to the underlying matrix */
			mpq_t* m;
			/** Row index of this iterator in the matrix */
			matrix_mpq::size_type i;
			/** Dimension of the matrix rows. */
			matrix_mpq::size_type d;
		};
		
		/** Constructor.
		 *  @param n			the number of rows in the matrix
		 *  @param d			the number of columns in the matrix
		 */
		matrix_mpq(size_type n, size_type d);
		
		/** Copy constructor. 
		 *  @param that			matrix to copy
		 */
		matrix_mpq(matrix_mpq const& that);
		
		/** Destructor */
		~matrix_mpq();
		
		/** Assignment operator.
		 *  @param that			matrix to copy into this one
		 */
		matrix_mpq& operator= (matrix_mpq const& that);
		
		/** Gets an iterator at the first row of the matrix. Note that the 
		 *  vector_mpq returned by this iterator shares storage with this 
		 *  matrix; if a copy of a matrix row is required, it must be 
		 *  explicitly made.
		 */
		iterator begin();
		
		/** Gets a constant iterator at the first row of the matrix. Note that 
		 *  the vector_mpq returned by this iterator shares storage with this 
		 *  matrix; if a copy of a matrix row is required, it must be 
		 *  explicitly made. */
		const_iterator begin() const;
		
		/** Gets an iterator one past the last row of the matrix. Note that the 
		 *  vector_mpq returned by this iterator shares storage with this 
		 *  matrix; if a copy of a matrix row is required, it must be 
		 *  explicitly made. */
		iterator end();
		
		/** Gets a constant iterator one past the last row of the matrix. Note 
		 *  that the vector_mpq returned by this iterator shares storage with 
		 *  this matrix; if a copy of a matrix row is required, it must be 
		 *  explicitly made. */
		const_iterator end() const;
		
		/** Gets the number of rows in this matrix. */
		size_type size() const;
		
		/** Gets the dimension of the rows in this matrix. */
		size_type dim() const;
		
		/** Indexing operator. Returns a mutable row reference. Note that the 
		 *  vector_mpq returned by this operator shares storage with this 
		 *  matrix; if a copy of a matrix row is required, it must be 
		 *  explicitly made.
		 *  @param i		The index of the row to return.
		 */
		reference operator[] (size_type i);
		
		/** Indexing operator. Returns a constant row reference. Note that the 
		 *  vector_mpq returned by this operator shares storage with this 
		 *  matrix; if a copy of a matrix row is required, it must be 
		 *  explicitly made.
		 *  @param i		The index of the row to return.
		 */
		const_reference operator[] (size_type i) const;
		
		/** Prints the matrix with space-separated rows inside square brackets.
		 *  @param o		The output stream to print on
		 *  @param m		The matrix to print
		 *  @return the same output stream
		 */
		friend std::ostream& operator<< (std::ostream& o, matrix_mpq const& m);
		
		/** Less than. Equivalent to a.compare(b) < 0 */
		friend bool operator<  (matrix_mpq const& a, matrix_mpq const& b);
		/** Equal. Equivalent to a.compare(b) == 0 */
		friend bool operator== (matrix_mpq const& a, matrix_mpq const& b);
		/** Greater than. Equivalent to a.compare(b) > 0 */
		friend bool operator>  (matrix_mpq const& a, matrix_mpq const& b);
		/** Less than or equal. Equivalent to a.compare(b) <= 0 */
		friend bool operator<= (matrix_mpq const& a, matrix_mpq const& b);
		/** Not equal. Equivalent to a.compare(b) != 0 */
		friend bool operator!= (matrix_mpq const& a, matrix_mpq const& b);
		/** Greater than or equal. Equivalent to a.compare(b) >= 0 */
		friend bool operator>= (matrix_mpq const& a, matrix_mpq const& b);
	
	private:
		/** Compares this matrix to another, lexicographically by rows.
		 *  @param that		The matrix to compare to
		 *  @return -1 for this matrix less than that, 0 for this matrix 
		 *  		equal to that, or 1 for this matrix greater than that.
		 */
		int compare(matrix_mpq const& that) const;
		
		/** Matrix data storage */
		mpq_t* m;
		/** Number of rows in the matrix */
		size_type n;
		/** Dimension of the matrix rows */
		size_type d;
	};
	
} /* namespace lrs */

#endif /* _LRS_MATRIX_HPP_ */