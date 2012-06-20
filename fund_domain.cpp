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

	fund_domain::const_iterator fund_domain::begin() const { return p.begin(); }

	fund_domain::const_iterator fund_domain::end() const { return p.end(); }

} /* namespace basil */
