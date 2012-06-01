/** Implements fundamental domain calculations from fund_domain.hpp.
 *
 *  @author Aaron Moss
 */

#include <vector>

#include <gmpxx.h>

#include "basil.hpp"
#include "fund_domain.hpp"

namespace basil {

	fund_domain::fund_domain(matrix_mpq qInv) : qInv(qInv) {
		p = std::vector<vector_mpq>();
	}

	void fund_domain::add_constraint(fund_domain::vector_mpq const& a,
			fund_domain::vector_mpq const& b) {

		vector_mpq h = row_mat_mul((a - b), qInv);
		p.push_back(h);
	}

	bool fund_domain::contains(fund_domain::vector_mpq const& x) const {

		ind n = size(), d = dim();

		/* for each row p_i, checks that p_i[0] + p_i[1..] * x >= 0 */
		for (ind i = 0; i < n; ++i) {
			mpq_class a = p[i][0];
			for (ind j = 0; j < d; ++j) {
				a += p[i][j+1] * x[j];
			}
			if ( sgn(a) < 0 ) return false;
		}

		return true;
	}

	ind fund_domain::dim() const { return qInv.dim(); }

} /* namespace basil */
