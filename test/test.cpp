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
	G() : z() {
		lrs::matrix_mpq m(4,3);
		int m_i[4][3] = {{1,-1,0},{1,0,-1},{1,0,1},{1,1,0}};
		for (long i = 0; i < 4; ++i) for (long j = 0; j < 3; ++j) {
			m[i][j] = mpq_class(m_i[i][j]);
		}
// 		std::cout << "\nm:" << m;
		g = basil::constructGram(m);
	}
	
	~G() {}
	
	basil::gram_matrix g;
	basil::gram_matrix z;
};

/** @return string representation of gram matrix */
std::string str(basil::gram_matrix const& g) {
	std::ostringstream o;
	o << g;
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

BOOST_AUTO_TEST_SUITE_END() /* gram_suite */

struct P {
	P() : factorizer() {}
	
	~P() {}
	
	basil::factor_list factorize(int x) {
		return *factorizer( mpz_class(x) );
	}
	
	int defactorize(basil::factor_list const& l) {
		return factorizer( l ).get_si();
	}
	
	basil::prime_factorizer factorizer;
};

BOOST_FIXTURE_TEST_SUITE( prime_suite, P )

BOOST_AUTO_TEST_CASE( factorization_test ) {
	using basil::factor_list;
	
	/* factorizing 0 gives empty pointer */
	BOOST_CHECK( ! factorizer( mpz_class(0) ) );
	
	/* factorizing 1 gives an empty factor list */
	factor_list l1;
	BOOST_CHECK( factorize(1) == l1 );
	BOOST_CHECK( defactorize( l1 ) == 1 );
	
	int i3[] = {0, 1};
	factor_list l3(i3, i3+2);
	BOOST_CHECK( factorize(3) == l3 );
	BOOST_CHECK( defactorize( l3 ) == 3 );
	
	int i10[] = {1, 0, 1};
	factor_list l10(i10, i10+3);
	BOOST_CHECK( factorize(10) == l10 );
	BOOST_CHECK( defactorize( l10 ) == 10 );
	
	int i99[] = {0, 2, 0, 0, 1};
	factor_list l99(i99, i99+5);
	BOOST_CHECK( factorize(99) == l99 );
	BOOST_CHECK( defactorize( l99 ) == 99 );
	
	/* check again, to test that the cache works properly */
	BOOST_CHECK( factorize(99) == l99 );
}

BOOST_AUTO_TEST_SUITE_END() /* prime_suite */
