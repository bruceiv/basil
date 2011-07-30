#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE lrs_wrapper_test
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <sstream>
#include <string>

#include <gmpxx.h>

#include "../gram.cpp" /* yes, I know this is evil - deal with it, it's the 
						* only way to get at the implementation details I put 
						* in here */
#include "../lrs/matrix.hpp"

/** Fixture for testing gram matrix helpers. */
struct G {
	G() : g(4), z() {
		int g_i[4][4] = {{1,2,2,0},{2,1,0,2},{2,0,1,2},{0,2,2,1}};
		for (long i = 0; i < 4; ++i) for (long j = 0; j < 4; ++j) {
			g(i,j) = g_i[i][j];
		}
	}
	
	~G() {}
	
	basil::gram_matrix g;
	basil::gram_matrix z;
};

/** @return string representation of the argument */
template <typename T>
std::string str(T const& x) {
	std::ostringstream o;
	o << x;
	return o.str();
}

BOOST_FIXTURE_TEST_SUITE( gram_suite, G )

BOOST_AUTO_TEST_CASE( dim_test ) {
	BOOST_CHECK( z.dim() == 0 );
	BOOST_CHECK( g.dim() == 4 );
}

BOOST_AUTO_TEST_CASE( equality_test ) {
// 	std::cout << "\nz: `" << str(z) << "'";
	BOOST_CHECK( str(z) == "|" );
	
// 	std::cout << "\ng: `" << str(g) << "'" << std::endl;
	BOOST_CHECK( str(g) == "| 1 2 2 0 | 2 1 0 2 | 2 0 1 2 | 0 2 2 1 |" );
}

BOOST_AUTO_TEST_CASE( construct_test ) {
	lrs::matrix_mpq m(4,3);
	int m_i[4][3] = {{1,-1,0},{1,0,-1},{1,0,1},{1,1,0}};
	for (long i = 0; i < 4; ++i) for (long j = 0; j < 3; ++j) {
		m[i][j] = mpq_class(m_i[i][j]);
	}
// 	std::cout << "\nm:" << m;
	basil::gram_matrix c = basil::constructGram(m);
	
	BOOST_CHECK( c == g );
}

BOOST_AUTO_TEST_SUITE_END() /* gram_suite */

struct P {
	P() : factorizer() {
		l1 = basil::factor_list();
		int i3[] = {0, 1};
		l3 = basil::factor_list(i3, i3+2);
		int i10[] = {1, 0, 1};
		l10 = basil::factor_list(i10, i10+3);
		int i12[] = {2, 1};
		l12 = basil::factor_list(i12, i12+2);
		int i22[] = {1, 0, 0, 0, 1};
		l22 = basil::factor_list(i22, i22+5);
		int i30[] = {1, 1, 1};
		l30 = basil::factor_list(i30, i30+3);
		int i99[] = {0, 2, 0, 0, 1};
		l99 = basil::factor_list(i99, i99+5);
		int i2970[] = {1, 3, 1, 0, 1};
		l2970 = basil::factor_list(i2970, i2970+5);
	}
	
	~P() {}
	
	basil::factor_list factorize(int x) {
		return factorizer( mpz_class(x) );
	}
	
	int defactorize(basil::factor_list const& l) {
		return factorizer( l ).get_si();
	}
	
	basil::prime_factorizer factorizer;
	
	basil::factor_list l1;
	basil::factor_list l3;
	basil::factor_list l10;
	basil::factor_list l12;
	basil::factor_list l22;
	basil::factor_list l30;
	basil::factor_list l99;
	basil::factor_list l2970;
};

BOOST_FIXTURE_TEST_SUITE( prime_suite, P )

BOOST_AUTO_TEST_CASE( factorization_test ) {
	
	/* factorizing 1 gives an empty factor list */
	BOOST_CHECK( factorize(1) == l1 );
	BOOST_CHECK( defactorize( l1 ) == 1 );
	
	BOOST_CHECK( factorize(3) == l3 );
	BOOST_CHECK( defactorize( l3 ) == 3 );
	
	BOOST_CHECK( factorize(10) == l10 );
	BOOST_CHECK( defactorize( l10 ) == 10 );
	
	BOOST_CHECK( factorize(99) == l99 );
	BOOST_CHECK( defactorize( l99 ) == 99 );
	
	/* check again, to test that the cache works properly */
	BOOST_CHECK( factorize(1) == l1 );
	BOOST_CHECK( factorize(99) == l99 );
}

BOOST_AUTO_TEST_CASE( mult_test ) {
	using basil::factor_list;
	using basil::mult;
	
	factor_list t1 = l1;
	BOOST_CHECK( mult(t1, l1) == l1 );
	/* t1 == l1 here, but make sure */
	t1 = l1;
	BOOST_CHECK( mult(t1, l3) == l3 );
	
	factor_list t3 = l3;
	BOOST_CHECK( mult(t3, l10) == l30 );
	/* check that t3 == t30 afterward */
	BOOST_CHECK( t3 == l30 );
	t3 = l30;
	BOOST_CHECK( mult(t3, l99) == l2970 );
}

BOOST_AUTO_TEST_CASE( mpr_test ) {
	using basil::mpr;
	
	mpr m0;
	mpr m1(mpz_class(1), mpz_class(1), mpz_class(1));
	mpr m2o3(mpz_class(2), mpz_class(1), mpz_class(3));
	mpr m2o3b(mpz_class(6), mpz_class(1), mpz_class(9));
	mpr m2o3c; m2o3c.n = 6; m2o3c.d = 9;
	mpr m2o3d(mpz_class(1), mpz_class(4), mpz_class(3));
	
	BOOST_CHECK( str(m0) == "0" );
	BOOST_CHECK( str(m1) == "1" );
	BOOST_CHECK( str(m2o3) == "2/3" );
	BOOST_CHECK( str(m2o3b) == "2/3" );
	BOOST_CHECK( str(m2o3c) == "6/9" );
	BOOST_CHECK( str(m2o3d) == "1r4/3" );
	
	BOOST_CHECK( m2o3 == m2o3b );
	BOOST_CHECK( m2o3 != m2o3c );
	BOOST_CHECK( m2o3 != m2o3d );
	
	m0 = m2o3;
	BOOST_CHECK( m0 == m2o3 );
	BOOST_CHECK( m0 == m2o3b );
}

BOOST_AUTO_TEST_CASE( norm_test ) {
	using basil::mpr;
	using basil::factor_list;
	
	mpq_class q(7,11);
	mpz_class ni(10);
	mpq_class nj(3);
	factor_list fi = l12; basil::mult( fi, l22 );
	factor_list fj = l12;
	
	mpr a = basil::norm(q, ni, nj, fi, fj, factorizer);
// std::cout << a << std::endl;
	mpr e(mpz_class(14), mpz_class(22), mpz_class(55));
	
	BOOST_CHECK( a == e );
}

BOOST_AUTO_TEST_SUITE_END() /* prime_suite */
