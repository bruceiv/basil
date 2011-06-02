/** Implements lrs::lrs class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdio>
#include <new>

#include "clrs.hpp"
#include "lrs.hpp"
#include "matrix.hpp"

namespace lrs {
	
	lrs::lrs(const matrix& m) throw(std::bad_alloc) {
		/* Initialize LRS */
		lrs_init_quiet(stdin, stdout);
		
		/* Init LRS global data */
		Q = lrs_alloc_dat((char*)"LRS globals");
		if (Q == 0) throw std::bad_alloc();
		
		/* Init LRS LP dictionary */
		ind n = m.n(), d = m.d();
		
		Q->m = n; Q->n = d;
		Q->geometric = true;
		
		P = lrs_alloc_dic(Q);
		if (P == 0) throw std::bad_alloc();
		
		matrix& m_nc = const_cast<matrix &>(m);
		for (ind i = 0; i < n; i++) {
			lrs_set_row_mp(P, Q, i+1, m_nc[i].num(), m_nc[i].den(), ge );
		}
	}
	
	lrs::~lrs() {
		/* FIXME This one I'm not quite sure about the memory management for */
		//if (Q->nredundcol > 0) lrs_clear_mp_matrix(Lin, Q->nredundcol, Q->n);
		lrs_free_dic(P, Q);
		lrs_free_dat(Q);
		lrs_close_quiet();
	}
	
	bool lrs::getFirstBasis() {
		/* Lin is an out parameter of this method, so it isn't 'initialized */
		return lrs_getfirstbasis(&P, Q, &Lin, true);
	}
	
	void lrs::printDict()
	{
		printA(P, Q);
	}
	
} /* namespace lrs */
