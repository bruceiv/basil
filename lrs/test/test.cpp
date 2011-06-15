#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE lrs_wrapper_test
#include <boost/test/unit_test.hpp>

#include "../cobasis.hpp"

/** Fixture for testing index_set helpers. */
struct C {
	C() : s1(16, 0xffe0ul ), s2(16, 0x007ful ) {}
	
	~C() {}
	
	lrs::index_set s1;
	lrs::index_set s2;
}

BOOST_FIXTURE_TEST_SUITE( lrs_wrapper_suite, C )

BOOST_AUTO_TEST_CASE( index_set_iter_test ) {
	
}

BOOST_AUTO_TEST_SUITE_END() /* lrs_wrapper_suite */
