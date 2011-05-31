/** Implements lrs::lrs class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdio>
#include <cstdlib>

#include "lrs.hpp"

#include "lrslib.h"


namespace lrs {
	
	lrs::lrs() {
		lrs_init_quiet(stdin, stdout);
		Q = lrs_alloc_dat((char*)"LRS globals");
	}
	
	lrs::~lrs() {
		free(Q);
		lrs_close_quiet();
	}
	
	bool lrs::loadMatrix(const matrix& m) {
		long n = m.n, d = m.d;
		
		Q->m = n; Q->n = d;
		P = lrs_alloc_dic(Q);
		if (P  == 0) return false;
		
		matrix& m_nc = const_cast<matrix &>(m);
		for (ind i = 0; i < n; i++) {
			lrs_set_row_mp(P, Q, i+1, m_nc[i].num(), m_nc[i].den(), ge );
		}
		
		return true;
	}

	
} /* namespace lrs */
