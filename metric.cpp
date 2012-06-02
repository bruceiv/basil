/** Implements distance metric computations from metric.hpp.
 *
 *  @author Aaron Moss
 */

#include <gmpxx.h>

#include "basil.hpp"
#include "metric.hpp"
#include "lrs/matrix.hpp"
#include "prime/factors.hpp"

namespace basil {

	////////////////////////////////////////////////////////////////////////////
	// mpr and matrix_mpr classes
	////////////////////////////////////////////////////////////////////////////

	mpr::mpr() : n(0), r(1), d(1) {}

	mpr::mpr(mpz_class n_, mpz_class r_, mpz_class d_) : n(n_), r(r_), d(d_) {
		/* reduce the numerator and denominator to lowest terms */
		mpz_class g;
		mpz_gcd(g.get_mpz_t(), n.get_mpz_t(), d.get_mpz_t());
		mpz_divexact(n.get_mpz_t(), n.get_mpz_t(), g.get_mpz_t());
		mpz_divexact(d.get_mpz_t(), d.get_mpz_t(), g.get_mpz_t());
	}

	bool operator== (mpr const& a, mpr const& b)
		{ return a.n == b.n && a.r == b.r && a.d == b.d; }

	bool operator!= (mpr const& a, mpr const& b)
		{ return !(a == b); }

	std::ostream& operator<< (std::ostream& o, mpr x) {
		o << x.n;
		if ( x.r != 1 ) o << "r" << x.r;
		if ( x.d != 1 ) o << "/" << x.d;

		return o;
	}

	matrix_mpr::matrix_mpr(ind n, ind d) : m(new mpr[n*d]), n(n), d(d) {}

	matrix_mpr::matrix_mpr(matrix_mpr const& that)
			: m(new mpr[that.n*that.d]), n(that.n), d(that.d) {
		for (ind i = 0; i < n*d; i++) m[i] = that.m[i];
	}

	matrix_mpr::~matrix_mpr() {
		delete[] m;
	}

	matrix_mpr& matrix_mpr::operator= (matrix_mpr const& that) {

		if (m != that.m) {
			if ( n*d != that.n*that.d ) {
				delete[] m;
				n = that.n;
				d = that.d;
				m = new mpq_class[n*d];
			}

			for (ind i = 0; i < n*d; i++) m[i] = that.m[i];
		}

		return *this;
	}

	ind matrix_mpr::size() const { return n; }

	ind matrix_mpr::dim() const { return d; }

	mpr& matrix_mpr::elem(ind i, ind j) { return m[i*d+j]; }

	mpr const& matrix_mpr::elem(ind i, ind j) const { return m[i*d+j]; }

	////////////////////////////////////////////////////////////////////////////
	// Metric matrices and helpers
	////////////////////////////////////////////////////////////////////////////

	lrs::matrix_mpq innerProdMat(lrs::matrix_mpq const& M) {

		ind n = M.size();

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

	lrs::matrix_mpq invQMat(lrs::matrix_mpq const& M) {

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

		return lu_inv(Q);
	}

	lrs::matrix_mpq orthoAugment(lrs::matrix_mpq const& M, bool augSigned) {

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

	lrs::matrix_mpq transformedInnerProdMat(lrs::matrix_mpq const& M,
									   lrs::matrix_mpq const& T) {
		ind n = M.size();

		lrs::matrix_mpq P(n,n);

		for (ind i = 0; i < n; ++i) {
			lrs::vector_mpq w = row_mat_mul(M.row(i), T);
			for (ind j = 0; j < n; ++j) P.elem(i, j) = inner_prod(w, M.row(j));
		}

		return P;
	}

	/** Calculates the inner product ip*sqrt(fi*fj)/(ni*nj), normalized
	 *  according to the standard rules for an mpr (numerator and denominator
	 *  share no factors, radical has no factors which are perfect squares).
	 */
	mpr norm(mpq_class ip, mpz_class ni, mpz_class nj, prime::factor_list fi,
			prime::factor_list fj, prime::factorizer product) {

		/* if ip == 0, return default (which initializes to 0, in canonical
		 * form) */
		if ( ip == 0 ) return mpr();

		/* first, multiply fj by fi */
		mult(fj, fi);

		/* now factor all the squares out into fi. fj, at the end should only
		 * have 0 or 1 of each prime, and fi should have half of what was
		 * removed  */

		uind k;
		/* ensure fi is big enough */
		for (k = fi.size(); k < fj.size(); ++k) fi.push_back(0);
		/* perform square root */
		for (k = 0; k < fj.size(); ++k) {
			if ( fj[k] & 0x1 /* fj[k] odd */ ) {
				fi[k] = (fj[k]-1)/2;  fj[k] = 1;
			} else /* fj[k] even */ {
				fi[k] = fj[k]/2;      fj[k] = 0;
			}
		}

		/* having canonicalized the radical, now generate the return value in
		 * canonical form */
		return mpr( mpz_class( abs( ip.get_num() ) * product(fi) ),
					product(fj),
					mpz_class( ip.get_den() * ni * nj ) );
	}

	matrix_mpr normedInnerProdMat(lrs::matrix_mpq const& M) {

		ind n = M.size();

		/* prime factorization functor */
		prime::factorizer factor;

		/* 1/||m[i]|| = sqrt(a_d*a_n)/a_n -- nums[i] = a_n, facs[i] = a_n*a_d */
		std::vector<mpz_class> nums(n);
		std::vector<prime::factor_list> facs(n);

		/* calculate norm information */
		mpq_class t;
		for (uind i = 0; i < n; ++i) {
			t = lrs::inner_prod(M.row(i), M.row(i));
			nums.push_back(t.get_num());
			/* NOTE: this assumes here that m[i] is not a zero vector - bad
			 * things happen otherwise */
			prime::factor_list fn = factor( t.get_num() );
			prime::factor_list fd = factor( t.get_den() );
			facs.push_back( mult(fn, fd) );
		}

		matrix_mpr P(n, n);

		/* calculate inner product matrix */
		for (uind i = 0; i < n; ++i) {
			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (uind j = 0; j < i; ++j) {

				t = lrs::inner_prod(M.row(i), M.row(j));
				mpr ip;

				if ( sgn(t) != 0 ) {
					ip = norm(t, nums[i], nums[j], facs[i], facs[j], factor);
				}

				P.elem(i, j) = ip;
			}
		}

		return P;
	}

} /* namespace basil */
