/** Implements distance metric computations from metric.hpp.
 *
 *  @author Aaron Moss
 */

#include <gmpxx.h>

#include "basil.hpp"
#include "metric.hpp"
#include "lrs/matrix.hpp"
#include "prime/factors.hpp"

#include "fmt.hpp"

namespace basil {

	////////////////////////////////////////////////////////////////////////////
	// mpr and matrix_mpr classes
	////////////////////////////////////////////////////////////////////////////

	/** Factors all the squares out of rf, placing their sqare roots in nf
	 *  @param nf				The numerator factors - will be set to the
	 *  						square root of the radical factors (should be
	 *  						at least as long as rf)
	 *  @param rf				The radical factors - will be reduced so that
	 *  						no perfect squares remain
	 */
	void sqrt(prime::factor_list& nf, prime::factor_list& rf) {
		for (uind k = 0; k < rf.size(); ++k) {
			if ( rf[k] & 0x1 /* rf[k] odd */ ) {
				nf[k] = rf[k]/2;	rf[k] = 1;
			} else /* rf[k] even */ {
				nf[k] = rf[k]/2;	rf[k] = 0;
			}
		}
	}

	mpr::mpr() : n(0), r(1), d(1) {}

	mpr::mpr(mpz_class n_, mpz_class r_, mpz_class d_) : n(n_), r(r_), d(d_) {}

	mpr mpr::makeNorm(mpz_class n_, prime::factor_list rf, mpz_class d_,
					  prime::factorizer factor) {

		/* factors to include in the numerator */
		prime::factor_list nf = prime::factor_list(rf.size());

		/* factor squares out of radical */
		sqrt(nf, rf);

		/* prepare return and normalize rational */
		mpr x(mpz_class(n_ * factor(nf)), factor(rf), d_);
		x.normRational();

		return x;
	}

	mpr::mpr(int x) : n(x), r(1), d(1) {}

	mpr& mpr::operator= (int x) { n = x; r = 1; d = 1; return *this; }

	void mpr::normRational() {
		/* get the GCD */
		mpz_class g;
		mpz_gcd(g.get_mpz_t(), n.get_mpz_t(), d.get_mpz_t());
		/* correct the sign */
		if ( sgn(d) < 0 ) mpz_neg(g.get_mpz_t(), g.get_mpz_t());
		/* factor out the GCD */
		if ( g != 1 ) {
			mpz_divexact(n.get_mpz_t(), n.get_mpz_t(), g.get_mpz_t());
			mpz_divexact(d.get_mpz_t(), d.get_mpz_t(), g.get_mpz_t());
		}
	}

	void mpr::norm(prime::factorizer factor) {
		/* factors of the radical */
		prime::factor_list rf = factor(r);
		/* factors to include in the numerator */
		prime::factor_list nf = prime::factor_list(rf.size());

		/* factor squares out of radical */
		sqrt(nf, rf);

		/* refactor numerator and radical */
		n *= factor(nf);
		r = factor(rf);
		/* re-normalize numerator and denominator */
		normRational();
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

	mpr abs(mpr const& x) { return mpr(abs(x.n), x.r, x.d); }

	int sgn(mpr const& x) { return sgn(x.n); }

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
				m = new mpr[n*d];
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

	lrs::matrix_mpq fixPlane(lrs::matrix_mpq const& M) {
		ind n = M.size(), d = M.dim();

		lrs::matrix_mpq F(n+1, d);

		for (ind i = 0; i < n; ++i) {
			F.row(i) = M.row(i);
		}

		lrs::vector_mpq v(d); v[0] = 1;
		F.row(n) = v;

		return F;
	}

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

//std::cout << "\tinnerProdMat():";
//for (ind i = 0; i < n; ++i) {
//std::cout << "\n\t";
//for (ind j = 0; j < n; ++j) std::cout << " " << P.elem(i, j);
//} std::cout << std::endl;

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

//std::cout << "\tQMat():";
//for (ind i = 0; i < d; ++i) {
//std::cout << "\n\t";
//for (ind j = 0; j < d; ++j) std::cout << " " << Q.elem(i, j);
//} std::cout << std::endl;

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

//std::cout << "\torthoAugment():";
//for (ind i = 0; i < n+rowAug; ++i) {
//std::cout << "\n\t";
//for (ind j = 0; j < d; ++j) std::cout << " " << R.elem(i, j);
//} std::cout << std::endl;

		return R;
	}

	lrs::matrix_mpq colRankAugment(lrs::matrix_mpq const& M,
								   lrs::index_set const& rows) {

		/* restrict matrix to selected rows */
		lrs::matrix_mpq B = M.row_restriction(rows);

		/* Find linearly independent columns of B */
		index_set colBasis = trans(B).lin_indep_rows();
		index_set missingCols = ~colBasis;
		missingCols.set(0, false);

		/* Augment the missing columns */
		lrs::matrix_mpq R(B.size() + missingCols.count(), B.dim());

		ind i = 0;
		for (; i < B.size(); ++i) {
			R.row(i) = B.row(i);
		}
		for (lrs::index_set_iter iter = lrs::begin(missingCols);
				iter != lrs::end(missingCols); ++iter) {
			R.elem(i, (*iter)-1) = 1;

			++i;
		}

//std::cout << "\t\tcolRankAugment():";
//for (ind i = 0; i < R.size(); ++i) {
//std::cout << "\n\t\t";
//for (ind j = 0; j < R.dim(); ++j) std::cout << " " << R.elem(i, j);
//} std::cout << std::endl;

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

//std::cout << "\ttransformedInnerProdMat():";
//for (ind i = 0; i < n; ++i) {
//std::cout << "\n\t";
//for (ind j = 0; j < n; ++j) std::cout << " " << P.elem(i, j);
//} std::cout << std::endl;

		return P;
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
		for (ind i = 0; i < n; ++i) {
			t = lrs::inner_prod(M.row(i), M.row(i));
			nums[i] = t.get_num();
			/* NOTE: this assumes here that m[i] is not a zero vector - bad
			 * things happen otherwise */
			prime::factor_list fn = factor( t.get_num() );
			prime::factor_list fd = factor( t.get_den() );
			facs[i] =  prime::mult(fn, fd);
		}

		matrix_mpr P(n, n);

		/* calculate inner product matrix */
		for (ind i = 0; i < n; ++i) {
			/* p[i][i] = 1, by def'n normed inner product */
			P.elem(i, i) = 1;

			/* Optimized here: p[i][j] = p[j][i], by def'n inner product */
			for (ind j = 0; j < i; ++j) {

				t = lrs::inner_prod(M.row(i), M.row(j));
				mpr ip;

				if ( sgn(t) != 0 ) {
					mpz_class num = t.get_num();
					prime::factor_list rad = facs[i]; prime::mult(rad, facs[j]);
					mpz_class den = t.get_den() * nums[i] * nums[j];

					ip = mpr::makeNorm(num, rad, den, factor);
				}

				P.elem(i, j) = ip;
				P.elem(j, i) = ip;
			}
		}

//std::cout << "\tnormedInnerProdMat():";
//for (ind i = 0; i < n; ++i) {
//std::cout << "\n\t";
//for (ind j = 0; j < n; ++j) std::cout << " " << P.elem(i, j);
//} std::cout << std::endl;

		return P;
	}

	lrs::matrix_mpq select_rows(lrs::matrix_mpq const& M, index_list const& l) {
		ind n = l.size(), d = M.dim();

		lrs::matrix_mpq R(n, d);

		for (ind i = 0; i < n; ++i) {
			R.row(i) = M.row(l[i]-1);
		}

		return R;
	}

} /* namespace basil */
