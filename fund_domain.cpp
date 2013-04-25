/** Implements fundamental domain calculations from fund_domain.hpp.
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

#include <set>
#include <vector>

#include <gmpxx.h>

#include "basil.hpp"
#include "fund_domain.hpp"
#include "metric.hpp"
#include "perm_utils.hpp"

namespace basil {

	fund_domain::fund_domain() : qInv() {
		p = std::vector<vector_mpq>();
	}

	fund_domain::fund_domain(matrix_mpq qInv) : qInv(qInv) {
		p = std::vector<vector_mpq>();
	}

	void fund_domain::build_from_seed(fund_domain::vector_mpq const& s,
			index_set const& sBasis,
			fund_domain::matrix_mpq const& A, permutation_list const& l) {

		ind n = A.size(), d = A.dim() - 1;

		/* augment input basis and matrix with x_0 = 1 plane (considered in
		 * homogenized coordinates). */
		index_list seedBasis = as_list(sBasis);
		index_list rowBasis(seedBasis);
		rowBasis.push_back(n+1);
		matrix_mpq Aa = fixPlane(A);

		/* Get first part of basis change matrix to transform out existing
		 * basis */
		matrix_mpq B = inv(select_rows(Aa, rowBasis));

		/* Find the images of the seed vertex under each permutation */
		std::set<vector_mpq> images;
		for (permutation_list::const_iterator iter = l.begin();
				iter != l.end(); ++iter) {

			permutation& p = **iter;

			/* calculate other part of basis change transformation matrix for
			 * permutation */
			index_list pBasis = apply(p, seedBasis);
			pBasis.push_back(n+1);
			matrix_mpq T = B * select_rows(Aa, pBasis);

			/** Get transformed vertex, ignoring permutations that fix the seed
			 *  vertex and vertices that have been seen before */
			vector_mpq v = mat_col_mul(T, s);
			if ( v != s ) {
				std::pair<std::set<vector_mpq>::iterator, bool> isNew
					= images.insert(v);
				if ( isNew.second ) {
					vector_mpq c = get_constraint(s, v);
std::cout << "Constraint " << c << " for image " << v << " of permutation " << p << std::endl;
					push_back(c);
				}
			}
		}

	}

	fund_domain::vector_mpq leading_unit(fund_domain::vector_mpq const& v) {
		ind n = v.size();

		fund_domain::vector_mpq r(n);
		ind i = 0;
		while ( i < n && v[i] == 0 ) {
			r[i] = 0;
			++i;
		}
		if ( i < n ) {
			mpq_class scale = abs(v[i]);
			while ( i < n ) {
				r[i] = v[i] / scale;
				++i;
			}
		}
		return r;
	}

	fund_domain::vector_mpq fund_domain::get_constraint(
			fund_domain::vector_mpq const& a,
			fund_domain::vector_mpq const& b) const {

		return leading_unit(row_mat_mul((a - b), qInv));
	}

	void fund_domain::push_back(fund_domain::vector_mpq const& c) {
		p.push_back(c);
	}

	bool fund_domain::contains(fund_domain::vector_mpq const& x) const {

		ind n = size(), d = dim();

		/* for each row p_i, checks that p_i * x >= 0 */
		for (ind i = 0; i < n; ++i) {
			mpq_class a = 0;
			for (ind j = 0; j < d; ++j) {
				a += p[i][j] * x[j];
			}
			if ( sgn(a) < 0 ) {
std::cout << "Vertex " << x << " eliminated by constraint " << p[i] << std::endl;
				return false;
			}
		}

		return true;
	}

	uind fund_domain::dim() const { return qInv.dim(); }

	uind fund_domain::size() const { return p.size(); }

	std::vector<fund_domain::vector_mpq> const&
			fund_domain::constraints() const { return p; }

	fund_domain::const_iterator fund_domain::begin() const { return p.begin(); }

	fund_domain::const_iterator fund_domain::end() const { return p.end(); }

} /* namespace basil */
