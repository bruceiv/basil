/** Implements lrs::lrs class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdio>
#include <new>
#include <sstream>

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
		/* Ported from David Bremner's equivalent code in lrsserv.c */
		
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
				if ( (col == 0L) || zero( A[Row[i]][col] ) ) {
					extraIncidences[incCount] = inequality[B[i] - lastdv];
					nIncidence++;
					incCount++;
				}
			}
		}
		
		mpz_class det(P->det);
		index_list cobInd(tempArray, tempArray + d);
		index_list extraInc(extraIncidences, extraIncidences + incCount);
		
		cobasis* cob = new cobasis(det, rflag, cobInd, nIncidence, extraInc);
		delete[] extraIncidences;
		return cob;
	}
	
	bool lrs::getFirstBasis() {
		/* Lin is an out parameter of this method, so it isn't initialized */
		return lrs_getfirstbasis(&P, Q, &Lin, true);
	}
	
	ind lrs::getRealDim() {
		return P->d;
	}
	
	vector_mpz* lrs::getSolution(ind col) {
		
		if (col < 0 || col > P->d) {
			std::ostringstream err;
			err << "getSolution: illegal column " << col;
			throw lrs_error( err.str() );
		}
		
		vector_mpz* output_p = new vector_mpz(Q->n);
		
		if (! lrs_getsolution(P, Q, output_p->v, col) ) {
			delete output_p;
			output_p = 0;
		}
		
		return output_p;
	}
	
	vector_mpz* lrs::getVertex() {
		ind i  = 1, iRedund = 0;
		
		ind nRedundCol = Q->nredundcol;
		ind* redundCol = Q->redundcol;
		ind n = Q->n;
		
		vector_mpz* output_p = new vector_mpz(n);
		vector_mpz& output = *output_p;
		
		/* copy column 0 to output */
		copy( output[0], P->det );
		
		for (ind j = 0; j < n; j++) {
			if (iRedund < nRedundCol && redundCol[iRedund] == j) {
				/* column was deleted as redundant */
				itomp(0L, output[j]);
				iRedund++;
			} else {
				/* column not deleted as redundant */
				getnextoutput(P, Q, i, 0L, output[j]);
				i++;
			}
		}
		
		reducearray(output.v, n);
		
		return output_p;
	}
	
	void lrs::printDict() {
		printA(P, Q);
	}
	
	void lrs::setCobasis(index_list& cob) {
		ind nlinearity = Q->nlinearity;
		ind* linearity = Q->linearity;
		ind* facet = Q->facet;
		ind m = Q->m;
		ind d = P->d;
		
		for (ind j = nlinearity, k = 0; j < d; j++, k++) {
			facet[j] = cob[k];
			
			/* check errors */
			if ( facet[j] < 1 || facet[j] > m ) {
				std::ostringstream err;
				err << "Start/restart cobasic indices must be in range [1," 
					<< m << "]";
				throw lrs_error(err.str());
			}
			for (ind i = 0; i < nlinearity; i++) {
				if ( facet[j] == linearity[i] ) {
					throw lrs_error(
						"Start/restart cobasic indices should not include "
						"linearities");
				}
			}
			for (ind i = 0; i < j; i++) {
				if ( facet[i] == facet[j] ) {
					
					throw lrs_error(
						"Start/restart cobasic indices must be distinct");
				}
			}
		}
			
		Q->restart=true;
		if ( !restartpivots(P, Q) )
			throw lrs_error("Could not restart pivots from given cobasis");
		
	}
	
} /* namespace lrs */
