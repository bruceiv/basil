/** Implmentation of "Gram" matrix computations from gram.hpp.
 * 
 *  @author Aaron Moss
 */

/*  Copyright: Aaron Moss, 2012, moss.aaron@unb.ca  */

/*  This file is part of Basil.

    Basil is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Basil is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Basil.  If not, see <http://www.gnu.org/licenses/>.  */

#include <algorithm>
#include <cstdlib>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include <boost/functional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include "basil.hpp"
#include "gram.hpp"
#include "metric.hpp"

namespace basil {
	
	gram_matrix::gram_matrix(uind n, uind k_) : n(n), k_(k_), m(new int*[n]), 
			m_(new int[n*n]) { 
		for (uind i = 0; i < n; ++i) m[i] = m_+(i*n);
		std::fill(m_, m_+(n*n), 0); /* zerofill matrix */
	}
	
	gram_matrix::gram_matrix(gram_matrix const& that) 
			: n(that.n), k_(that.k_), m(new int*[that.n]), 
			m_(new int[that.n*that.n]) { 
		for (uind i = 0; i < n; ++i) /* copy row indices */
			m[i] = m_+(that.m[i]-that.m_);
		std::copy(that.m_, that.m_+(n*n), m_); /* copy matrix */
	}
		
	gram_matrix::~gram_matrix() {
		delete[] m_;
		delete[] m;
	}
	
	gram_matrix& gram_matrix::operator= (gram_matrix const& that) {
		if (m_ != that.m_) {
			if (n != that.n) {
				n = that.n;
				delete[] m_;
				delete[] m;
				m = new int*[n];
				m_ = new int[n*n];
			}
			
			for (uind i = 0; i < n; ++i) /* copy row indices */
				m[i] = m_+(that.m[i]-that.m_);
			std::copy(that.m_, that.m_+(n*n), m_); /* copy matrix */
			k_ = that.k_;
		}
		
		return *this;
	}
	
	int& gram_matrix::at(uind i, uind j) { return m[i][j]; }
	
	int gram_matrix::at(uind i, uind j) const { return m[i][j]; }
	
	uind gram_matrix::dim() const { return n; }
	
	gram_matrix gram_matrix::restriction(index_set s) const {
		gram_matrix r(s.count(), k_);
		
		uind i = 0;
		for (lrs::index_set_iter iterI = lrs::begin(s); 
				iterI != lrs::end(s); ++iterI) {
			uind j = 0;
			for (lrs::index_set_iter iterJ = lrs::begin(s); 
					iterJ != lrs::end(s); ++iterJ) {
				/* correct for 1-indexed index_set_iter */
				r.m[i][j] = m[(*iterI)-1][(*iterJ)-1];
				++j;
			}
			++i;
		}
		
		return r;
	}
	
	gram_matrix gram_matrix::abs() const {
		gram_matrix a(n);
		for (uind i = 0; i < n; ++i) for (uind j = 0; j < n; ++j) 
			a.m[i][j] = std::abs(m[i][j]);
		return a;
	}
	
	gram_matrix gram_matrix::doubled() const {
		gram_matrix d(2*n);
		for (uind i = 0; i < n; ++i) for (uind j = 0; j < n; ++j) {
			uind iPos = 2*i, jPos = 2*j; uind iNeg = iPos+1, jNeg = jPos+1;
			int xPos = m[i][j]; int xNeg = -xPos;
			d.m[iPos][jPos] = d.m[iNeg][jNeg] = xPos;
			d.m[iPos][jNeg] = d.m[iNeg][jPos] = xNeg;
		}
		return d;
	}
	
	gram_matrix gram_matrix::permlibCanon() const {
		
		typedef boost::unordered_map<int,int> rep_map;
		
		gram_matrix c(n);
		rep_map reps;
		
		for (uind i = 0; i < n; ++i) for (uind j = 0; j < n; ++j) {
			int val = m[i][j];
			rep_map::iterator iter = reps.find(val);
			int rep;
			if ( iter == reps.end() ) {
				rep = c.k_++;
				reps.insert(std::make_pair(val, rep));
			} else {
				rep = iter->second;
			}
			c.m[i][j] = rep;
		}
		
		return c;
	}
	
	/** Functional to lexicograpically compare two rows of a gram matrix */
	class gram_matrix_row_comparator 
			: public std::binary_function<int*, int*, bool> {
	public:
		gram_matrix_row_comparator(uind n = 0) : n(n) {}
		
		bool operator() (int* a, int* b) const 
			{ return std::lexicographical_compare(a, a+n, b, b+n); }
	
	private:
		/* length of the rows to compare */
		uind n;
	}; /* class gram_matrix_row_comparator */
	
	gram_matrix& gram_matrix::sort() {
		
		/* sort rows */
		for (uind i = 0; i < n; ++i) std::sort(m[i], m[i]+n);
		
		/* lex-sort matrix by rows */
		std::sort(m, m+n, gram_matrix_row_comparator(n));
		
		return *this;
	}
	
	bool operator== (gram_matrix const& a, gram_matrix const& b) { 
		bool isEqual = a.n == b.n;
		for (uind i = 0; isEqual && i < a.n; ++i) {
			isEqual = std::equal(a.m[i], a.m[i]+a.n, b.m[i]);
		}
		return isEqual;
	}
	
	std::ostream& operator<< (std::ostream& o, gram_matrix const& g) {
		for (uind i = 0; i < g.n; ++i) {
			o << "| ";
			for (uind j = 0; j < g.n; ++j) o << g.m[i][j] << " ";
		}
		o << "|";
		
		return o;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Hasher implementation
	////////////////////////////////////////////////////////////////////////////
	
	std::size_t gram_matrix_hash::operator() (gram_matrix const& m) const {
		std::size_t seed = 0UL;
		for (uind i = 0; i < m.n; ++i) for (uind j = 0; j < m.n; ++j) {
			/* combine values into hash */
			boost::hash_combine(seed, m.m[i][j] );
		}
		return seed;
	}
	
	/** Functional to hash an mpr. */
	class mpr_hash : public std::unary_function<mpr, std::size_t> {
	public:
		std::size_t operator() (mpr const& x) const {
			std::size_t seed = 0UL;
			/* combine low-order bits of values into hash */
			boost::hash_combine(seed, x.n.get_si());
			boost::hash_combine(seed, x.r.get_si());
			boost::hash_combine(seed, x.d.get_si());
			return seed;
		}
	}; /* class mpr_hash */
	
	/** Functional to hash an mpq_class */
	class mpq_class_hash : public std::unary_function<mpq_class, std::size_t> {
	public:
		std::size_t operator() (mpq_class const& x) const {
			std::size_t seed = 0UL;
			boost::hash_combine(seed, x.get_num().get_si());
			boost::hash_combine(seed, x.get_den().get_si());
			return seed;
		}
	}; /* class mpq_class_hash */

	/** Constructs a Gram matrix according to the rules for constructGram.
	 *  @param Elem			The element type. Needs a comparisons for equality,
	 *  					as well as abs(x) and sgn(x) defined for Elem x.
	 *  					Should also be assignable from int.
	 *  @param Mat			The matrix type. Needs methods size() and dim() to
	 *  					determine the number of rows and columns,
	 *  					respectively, as well as an accessor elem(i, j) to
	 *  					return the element (of type Elem) at row i, column j
	 *  @param Hash			A hashing functor of signature size_t Hash(Elem x)
	 */
	template<typename Elem, typename Mat, typename Hash>
	gram_matrix constructGram(Mat const& M) {
		/* typedefs */
		typedef typename boost::unordered_map<Elem, int, Hash> Elem_Map;
		typedef typename Elem_Map::iterator Elem_Iter;

		/* size of the matrix */
		ind n = M.size();

		/* get unique representatives */
		Elem_Map reps;
		Elem zero(0); reps.insert(std::make_pair(zero, 0));
		int nextRep = 1;

		gram_matrix G(n, 1);

		for (ind i = 0; i < n; ++i) for (ind j = 0; j < n; ++j) {
			int rep;
			Elem val = M.elem(i, j);
			Elem key = abs(val);

			Elem_Iter res = reps.find(key);
			if ( res == reps.end() ) {
				/* representative not found, make new */
				rep = nextRep++;
				reps.insert(std::make_pair(key, rep));
			} else {
				/* representative found, use */
				rep = res->second;
			}
			rep *= sgn(val);
			G(i, j) = rep;
		}

		return G;
	}

	gram_matrix constructGram(lrs::matrix_mpq const& M) {
		return constructGram<mpq_class, lrs::matrix_mpq, mpq_class_hash>(M);
	}

	gram_matrix constructGram(matrix_mpr const& M) {
		return constructGram<mpr, matrix_mpr, mpr_hash>(M);
	}

} /* namespace basil */
