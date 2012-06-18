/** Implements fundamental domain calculations from fund_domain.hpp.
 *
 *  @author Aaron Moss
 */

#include <vector>

#include <gmpxx.h>

#include "basil.hpp"
#include "fund_domain.hpp"

namespace basil {

	fund_domain::fund_domain() : qInv() {
		p = std::vector<vector_mpq>();
	}

	fund_domain::fund_domain(matrix_mpq qInv) : qInv(qInv) {
		p = std::vector<vector_mpq>();
	}

	void fund_domain::add_constraint(fund_domain::vector_mpq const& a,
			fund_domain::vector_mpq const& b) {

		vector_mpq h = row_mat_mul((a - b), qInv);
		p.push_back(a+h);
	}

	bool fund_domain::contains(fund_domain::vector_mpq const& x) const {

		ind n = size(), d = dim();

		/* for each row p_i, checks that p_i * x >= 0 */
		for (ind i = 0; i < n; ++i) {
			//mpq_class a = p[i][0];
			mpq_class a = 0;
			for (ind j = 0; j < d; ++j) {
				a += p[i][j] * x[j];
			}
			if ( sgn(a) < 0 ) return false;
		}

		return true;
	}

	uind fund_domain::dim() const { return qInv.dim(); }

	uind fund_domain::size() const { return p.size(); }

	std::vector<fund_domain::vector_mpq> const&
			fund_domain::constraints() const { return p; }

} /* namespace basil */
