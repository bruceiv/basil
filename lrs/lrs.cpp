/** Implements lrs::lrs class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdio>
#include <cstdlib>
#include <new>

#include "lrs.hpp"

#include "lrslib.h"


namespace lrs {
	
	lrs::lrs() throw(std::bad_alloc) {
		lrs_init_quiet(stdin, stdout);
		Q = lrs_alloc_dat((char*)"LRS globals");
		if (Q == 0) throw std::bad_alloc();
	}
	
	lrs::~lrs() {
		free(Q);
		lrs_close_quiet();
	}                                                                                                                 
	
	void lrs::loadMatrix(const matrix& m) throw(std::bad_alloc) {
		ind n = m.n(), d = m.d();
		
		Q->m = n; Q->n = d;
		Q->geometric = true;
		
		P = lrs_alloc_dic(Q);
		if (P  == 0) throw std::bad_alloc();
		
		matrix& m_nc = const_cast<matrix &>(m);
		for (ind i = 0; i < n; i++) {
			lrs_set_row_mp(P, Q, i+1, m_nc[i].num(), m_nc[i].den(), ge );
		}
	}
	
	void lrs::getFirstBasis() {
		/* FIXME - find out where/how to initialize linearities param */
		lrs_getfirstbasis(&P, Q, 0, true);
	}

	
} /* namespace lrs */
