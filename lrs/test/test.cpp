#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE lrs_wrapper_test
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>

#include <gmpxx.h>

#include "../cobasis.hpp"
#include "../matrix.hpp"

/** Fixture for testing index_set helpers. */
struct M {
	M() : m1(3, 4), m2(3, 4) {
		for (long i = 0; i < 3; ++i) for (long j = 0; j < 4; ++j) {
			m1[i][j] = mpq_class(10*i+j);
			m2.elem(i,j) = mpq_class(10*i+j);
		}
	}
	
	~M() {}
	
	lrs::matrix_mpq m1;
	lrs::matrix_mpq m2;
};

BOOST_FIXTURE_TEST_SUITE( matrix_suite, M )

BOOST_AUTO_TEST_CASE( equality_test ) {
	BOOST_CHECK( m1 == m2 );
	m1[2][2] = mpq_class(0);
	BOOST_CHECK( m1 != m2 );
	BOOST_CHECK( m1 < m2 );
	m2.elem(2,2) = mpq_class(0);
	BOOST_CHECK( m1 == m2 );
}

BOOST_AUTO_TEST_CASE( access_test ) {
	lrs::vector_mpq v(4);
	v[0] = mpq_class(20);
	v[1] = mpq_class(21);
	v[2] = mpq_class(22);
	v[3] = mpq_class(23);
	
	BOOST_CHECK( v == m1[2] );
	BOOST_CHECK( v == m2.row(2) );
}

BOOST_AUTO_TEST_CASE( inner_prod_test ) {
	/* [ [ 14, 74, 134 ], [ 74, 534, 994 ], [ 134, 994, 1854 ] ] */
	int ipi[3][3] = {{14,74,134},{74,534,994},{134,994,1854}};
	lrs::matrix_mpq ipm(3,3);
	for (long i = 0; i < 3; ++i) for (long j = 0; j < 3; ++j) {
		ipm[i][j] = mpq_class(ipi[i][j]);
	}
	
	BOOST_CHECK( ipm == m1.inner_prod_mat() );
	BOOST_CHECK( ipm == m2.inner_prod_mat() );
}

BOOST_AUTO_TEST_CASE( restriction_test ) {
	lrs::index_set rs(4); rs.set(1).set(3);
	
	lrs::matrix_mpq rm(2,2);
	rm[0][0] = mpq_class(0);  rm[0][1] = mpq_class(2);
	rm[1][0] = mpq_class(20); rm[1][1] = mpq_class(22);
	
	BOOST_CHECK( rm == m1.restriction(rs) );
	BOOST_CHECK( rm == m2.restriction(rs) );
}

BOOST_AUTO_TEST_CASE( sort_test ) {
	int v1i[4] = { 13, 11, 10, 12 };
	lrs::vector_mpq v1m(4);
	for (long i = 0; i < 4; ++i) v1m[i] = v1i[i];
	
	std::sort(v1m.begin(), v1m.end());
	
	BOOST_CHECK( v1m == m1[1] );
	BOOST_CHECK( v1m == m2.row(1) );
	
	int mui[3][4] = {{10,11,12,13},{20,21,22,23},{0,1,2,3}};
	lrs::matrix_mpq mum(3,4);
	for (long i = 0; i < 3; ++i) for (long j = 0; j < 4; ++j) 
		mum[i][j] = mui[i][j];
	
	std::sort(mum.begin(), mum.end());
	
	BOOST_CHECK( mum == m1 );
	BOOST_CHECK( mum == m2 );
}

BOOST_AUTO_TEST_CASE( assign_test ) {
	lrs::vector_mpq v1 = m1[1];
	lrs::vector_mpq v2 = m2.row(1);
	lrs::vector_mpq v3(m1[1]);
	
	BOOST_CHECK( v1 == m1[1] );
	BOOST_CHECK( v2 == m1.row(1) );
	BOOST_CHECK( v3 == v1 );
}

BOOST_AUTO_TEST_CASE( matrix_rearrange ) {
	lrs::matrix_mpq m3(3,4);
	
	m3[0] = m1[2];
	m3.row(1) = m1[1];
	m3[2] = m2.row(0);
	
	lrs::matrix_mpq m4(m3);
	
	lrs::vector_mpq tmp;
	tmp = m3[0];
	m3[0] = m3[2];
	m3[2] = tmp;
	
	BOOST_CHECK( m3 == m1 );
	
	lrs::matrix_mpq::iterator iter0 = m4.begin();
	lrs::matrix_mpq::iterator iter2 = iter0; iter2 += 2;
	
	//lrs::vector_mpq tmp2; tmp2 = *iter0;
	lrs::vector_mpq tmp2 = *iter0;	
	*iter0 = *iter2;
	*iter2 = tmp2;
	
	BOOST_CHECK( m4 == m3 );
	BOOST_CHECK( m4 == m1 );
}

BOOST_AUTO_TEST_CASE( matrix_iter ) {
	long i = 0;
	for (lrs::matrix_mpq::const_iterator iter = m1.begin(); 
			iter != m1.end(); ++iter, ++i) {
		BOOST_CHECK( *iter == m2[i] );
	}
}

BOOST_AUTO_TEST_CASE( matrix_self_sort ) {
	int mui[3][4] = {{12,10,13,11},{23,21,22,20},{0,2,1,3}};
	lrs::matrix_mpq mum(3,4);
	for (long i = 0; i < 3; ++i) for (long j = 0; j < 4; ++j) 
		mum[i][j] = mui[i][j];
	
	for (lrs::matrix_mpq::iterator iter = mum.begin(); 
			iter != mum.end(); ++iter) {
		std::sort( (*iter).begin(), (*iter).end() );
	}
	
	int mri[3][4] = {{10,11,12,13},{20,21,22,23},{0,1,2,3}};
	lrs::matrix_mpq mrm(3,4);
	for (long i = 0; i < 3; ++i) for (long j = 0; j < 4; ++j) 
		mrm[i][j] = mri[i][j];
	
	BOOST_CHECK( mum == mrm );
	
	/* sorts the matrix */
	
// 	std::cout << "initial state:\t" << mum << std::endl;
	
	lrs::matrix_mpq::iterator iter = mum.begin();
	lrs::matrix_mpq::iterator iter2 = iter;
	iter += 2;
	
	lrs::vector_mpq tmp = *iter2;
// 	std::cout << "\t`tmp = *iter2;` tmp == " << tmp << " *iter == " << *iter 
// 			<< " *iter2 == " << *iter2 << std::endl;
	
	*iter2 = *iter;
// 	std::cout << "\t`*iter2 = *iter;` tmp == " << tmp << " *iter == " << *iter 
// 			<< " *iter2 == " << *iter2 << std::endl;
	
	*iter = tmp;
// 	std::cout << "\t`*iter = tmp;` tmp == " << tmp << " *iter == " << *iter 
// 			<< " *iter2 == " << *iter2 << std::endl;
// 	
// 	std::cout << "first swap:\t" << mum << std::endl;
	
	++iter2;
	tmp = *iter2; *iter2 = *iter; *iter = tmp;
	
// 	std::cout << "final swap:\t" << mum << std::endl;
	
	BOOST_CHECK( mum == m1 );
}

BOOST_AUTO_TEST_SUITE_END() /* lrs_wrapper_suite */
