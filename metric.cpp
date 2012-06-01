#include "basil.hpp"
#include "metric.hpp"

#include "lrs/matrix.hpp"

namespace basil {

	lrs::matrix_mpq inner_prod_mat(lrs::matrix_mpq const& M) {

		ind n = M.size(), d = M.dim();

		lrs::matrix_mpq P(n,n);

		mpq_class t;
		for (ind i = 0; i < n; ++i) {
			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (ind j = 0; j < i; ++j) {
				t = inner_prod(M.row(i), M.row(j));
				P.elem(i,j) = t;
				P.elem(j,i) = t;
			}
			/* Handle j = i case here */
			t = inner_prod(M.row(i), M.row(i));
			P.elem(i,i) = t;
		}

		return P;
	}

	lrs::matrix_mpq q_mat(lrs::matrix_mpq const& M) {

		ind n = M.size(), d = M.dim();

		lrs::matrix_mpq Q(d, d);

		/* calculate upper triangle */
		for (ind i = 0; i < n; ++i) {
			for (ind j = 0; j < d; ++j) {
				for (ind k = j; k < d; ++k) {
					Q.elem(j, k) += M.elem(i, j) * M.elem(i, k);
				}
			}
		}

		/* copy into lower triangle */
		for (ind j = 1; j < d; ++j) for (ind k = 0; k < j; ++k) {
			Q.elem(j, k) = Q.elem(k, j);
		}

		return Q;
	}

	lrs::matrix_mpq ortho_augment(lrs::matrix_mpq const& M, bool augSigned) {

		using lrs::matrix_mpq;

		index_set goodRows = M.lin_indep_rows();

		ind n = M.size(), d = M.dim(), r = goodRows.count();

		matrix_mpq G = M.row_restriction(goodRows);

		index_set goodCols = trans(G).lin_indep_rows();
		index_set badCols = ~goodCols;
		badCols.set(0,false);

		matrix_mpq B = G.col_restriction(goodCols);
		matrix_mpq C = G.col_restriction(badCols);

		matrix_mpq A = inv(B) * -C;

		int rowAug = ( augSigned ) ? 2*(d-r) : (d-r);
		matrix_mpq R(n+rowAug, d);

		for (ind i = 0; i < n; ++i) {
			R.row(i) = M.row(i);
		}

		if ( augSigned ) {
			for (ind j = 0; j < (d-r); ++j) {
				for (ind i = 0; i < r; ++i) {
					mpq_class x = A.elem(i, j);
					R.elem(n+2*j, i) = x;
					R.elem(n+2*j+1, i) = -x;
				}
				//incorporate identity matrix augmentation
				R.elem(n+2*j, r+j) = 1;
				R.elem(n+2*j+1, r+j) = -1;
			}
		} else {
			for (ind j = 0; j < (d-r); ++j) {
				for (ind i = 0; i < r; ++i) {
					R.elem(n+j, i) = A.elem(i, j);
				}
				//incorporate identity matrix augmentation
				R.elem(n+j, r+j) = 1;
			}
		}

		return R;
	}

} /* namespace basil */
