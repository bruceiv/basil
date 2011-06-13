#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE lru_cache_test
#include <boost/test/unit_test.hpp>

#include "../cache.hpp"

struct F {
	F() : c(5) {
		c.insert(2);
		c.insert(4);
		c.insert(6);
	}
	~F() {}
	
	lru::cache<int> c;
};

BOOST_FIXTURE_TEST_SUITE( lru_cache_suite, F )

BOOST_AUTO_TEST_CASE( testLookup ) {
									/* [2 4 6] */
	BOOST_CHECK( c.lookup(2) );		/* [4 6 2] */
	BOOST_CHECK( c.lookup(4) );		/* [6 2 4] */
	BOOST_CHECK( c.lookup(6) );		/* [2 4 6] */
	BOOST_CHECK( ! c.lookup(5) );	/* [2 4 6] */
}

BOOST_AUTO_TEST_CASE( testInsert ) {
									/* [2 4 6] */
	BOOST_CHECK( ! c.lookup(7) );	/* [2 4 6] */
	BOOST_CHECK( ! c.insert(7) );	/* [2 4 6 7] */
	BOOST_CHECK( c.lookup(7) );		/* [2 4 6 7] */
}

BOOST_AUTO_TEST_CASE( testRemove ) {
									/* [2 4 6] */
	BOOST_CHECK( c.lookup(2) );		/* [4 6 2] */
	BOOST_CHECK( c.remove(2) );		/* [4 6] */
	BOOST_CHECK( ! c.lookup(2) );	/* [4 6] */
}

BOOST_AUTO_TEST_CASE( testOverflow ) {
									/* [2 4 6] */
	BOOST_CHECK( c.lookup(6) );		/* [2 4 6] */
	BOOST_CHECK( ! c.insert(3) );	/* [2 4 6 3] */
	BOOST_CHECK( ! c.insert(1) );	/* [2 4 6 3 1] */
	BOOST_CHECK( c.insert(4) );		/* [2 6 3 1 4] */
	BOOST_CHECK( c.lookup(2) );		/* [6 3 1 4 2] */
	BOOST_CHECK( ! c.insert(5) );	/* [3 1 4 2 5] 
									 * cache should overflow, dump 6 here */
	BOOST_CHECK( ! c.lookup(6) );	/* [3 1 4 2 5] */
	BOOST_CHECK( c.lookup(1) );		/* [3 4 2 5 1] */
	BOOST_CHECK( c.lookup(2) );		/* [3 4 5 1 2] */
	BOOST_CHECK( c.lookup(3) );		/* [4 5 1 2 3] */
	BOOST_CHECK( c.lookup(4) );		/* [4 1 2 3 4] */
	BOOST_CHECK( c.lookup(5) );		/* [1 2 3 4 5] */
}

BOOST_AUTO_TEST_SUITE_END() /* lru_cache_suite */
