/** Implements lrs::lrs class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdio>
#include <new>
#include <valarray>

#include "clrs.hpp"
#include "cobasis.hpp"
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
	
	cobasis* lrs::getCobasis(ind col) {
		using std::valarray;
		
		ind i;
		ind rflag = -1;  /* for finding inequality number of ray column */
		
		/* set local variables for structures */
		matrix_t& A = P->A;
		ind* B = P->B;
		ind* C = P->C;
		ind* Col = P->Col;
		ind* Row = P->Row;
		ind* inequality = Q->inequality;
		ind* tempArray = Q->temparray;
		ind d = P->d;
		ind lastdv = Q->lastdv;
		ind m = P->m;
		
		ind nIncidence;  /* count number of tight inequalities */
		ind incCount = 0;
		ind* extraIncidences = new ind[m];
		
		for (i = 0; i < d; i++) {
			if (Col[i] == col) rflag = tempArray[i]; /* look for ray index */
			
			tempArray[i] = inequality[C[i] - lastdv];
		}
		
		nIncidence = (col == 0) ? d : d-1;
		
		for (i = lastdv+1; i <= m; i++) {
			if ( zero( A[Row[i]][0] ) ) {
				if ( (col == ZERO) || zero( A[Row[i]][col] ) ) {
					extraIncidences[incCount] = inequality[B[i] - lastdv];
					nIncidence++;
					incCount++;
				}
			}
		}
		
		cobasis* cob = new cobasis(mpz_class(P->det), rflag, 
								   std::valarray<ind>(tempArray, d), 
								   nIncidence, 
								   std::valarray<ind>(extraIncidences, incCount)
  								);
		delete[] extraIncidences;
		return cob;
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
