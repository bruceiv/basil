#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE lrs_wrapper_test
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <sstream>
#include <string>

#include <gmpxx.h>

#include "../gram.hpp"
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
