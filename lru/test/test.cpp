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
	int a246[] = {2, 4, 6};
	int a264[] = {2, 6, 4};
									/* [2 4 6] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a246, a246 + 3 );
	
	BOOST_CHECK( c.lookup(4) );		/* [2 6 4] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a264, a264 + 3 );
		
	BOOST_CHECK( ! c.lookup(5) );	/* [2 6 4] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a264, a264 + 3 );
}

BOOST_AUTO_TEST_CASE( testInsert ) {
	/* The collection equality for lookup has been tested by testLookup */
	int a2467[] = {2, 4, 6, 7};
	
									/* [2 4 6] */
	BOOST_CHECK( ! c.lookup(7) );	/* [2 4 6] */
	
	BOOST_CHECK( ! c.insert(7) );	/* [2 4 6 7] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a2467, a2467 + 4 );
	
	BOOST_CHECK( c.lookup(7) );		/* [2 4 6 7] */
}

BOOST_AUTO_TEST_CASE( testRemove ) {
	int a46[] = {4, 6};
									/* [2 4 6] */
	BOOST_CHECK( c.lookup(2) );		/* [4 6 2] */
	
	BOOST_CHECK( c.remove(2) );		/* [4 6] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a46, a46 + 2 );
	
	BOOST_CHECK( ! c.lookup(2) );	/* [4 6] */
}

BOOST_AUTO_TEST_CASE( testOverflow ) {
	int a2463[] = {2, 4, 6, 3};
	int a24631[] = {2, 4, 6, 3, 1};
	int a26314[] = {2, 6, 3, 1, 4};
	int a63142[] = {6, 3, 1, 4, 2};
	int a31425[] = {3, 1, 4, 2, 5};
	
									/* [2 4 6] */
	BOOST_CHECK( c.lookup(6) );		/* [2 4 6] */
	
	BOOST_CHECK( ! c.insert(3) );	/* [2 4 6 3] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a2463, a2463 + 4 );
	
	BOOST_CHECK( ! c.insert(1) );	/* [2 4 6 3 1] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a24631, a24631 + 5 );
	
	BOOST_CHECK( c.insert(4) );		/* [2 6 3 1 4] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a26314, a26314 + 5 );
	
	BOOST_CHECK( c.lookup(2) );		/* [6 3 1 4 2] */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a63142, a63142 + 5 );
	
	BOOST_CHECK( ! c.insert(5) );	/* [3 1 4 2 5] 
									 * cache should overflow, dump 6 here */
	BOOST_CHECK_EQUAL_COLLECTIONS( c.begin(), c.end(), a31425, a31425 + 5 );
	
	BOOST_CHECK( ! c.lookup(6) );	/* [3 1 4 2 5] */
}

BOOST_AUTO_TEST_SUITE_END() /* lru_cache_suite */
