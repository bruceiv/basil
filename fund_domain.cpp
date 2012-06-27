/** Implements fundamental domain calculations from fund_domain.hpp.
 *
 *  @author Aaron Moss
 */

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

		/* pre-processing for transformation matrices for permutations */
		//lrs::index_set rowBasis = A.lin_indep_rows();
		matrix_mpq B = inv(colRankAugment(A, sBasis));

		/* Find the images of the seed vertex under each permutation */
		std::set<vector_mpq> images;
		for (permutation_list::const_iterator iter = l.begin();
				iter != l.end(); ++iter) {

			permutation& p = **iter;

			/* get transformation matrix for permutation */
			lrs::index_set permBasis = apply(p, sBasis);
			matrix_mpq T = B * colRankAugment(A, permBasis);

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

		/* Add the constraint for each image */
//		for (std::set<vector_mpq>::iterator jter = images.begin();
//				jter != images.end(); ++jter) {
//			push_back(get_constraint(s, *jter));
//			vector_mpq c = get_constraint(s, *jter);
//			mpq_class a = 0;
//			for (ind j = 0; j < c.size(); ++j) {
//				a += c[j] * s[j];
//			}
//			if ( sgn(a) < 0 ) { std::cout << "WARNING: constraint " << c << " excludes seed " << s << std::endl; c = -c; }
//			push_back(c);
//		}
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
