/** Implements lrs::lrs class from lrs.hpp C++ wrapper for LRS.
 *
 *  @author Aaron Moss
 */

#include <cstdio>
#include <new>
#include <sstream>

#include <gmpxx.h>

#include "clrs.hpp"
#include "cobasis.hpp"
#include "lrs.hpp"
#include "matrix.hpp"

namespace lrs {
	
	lrs::lrs(matrix const& m, lrs_opts o) throw(std::bad_alloc) : o(o) {
		/* Initialize LRS */
		lrs_init_quiet(stdin, stdout);
		
		/* Init LRS global data */
		Q = lrs_alloc_dat((char*)"LRS globals");
		if (Q == 0) throw std::bad_alloc();
		
		/* Init LRS LP dictionary */
		initDat(Q, m.n(), m.d());
		
		P = lrs_alloc_dic(Q);
		if (P == 0) throw std::bad_alloc();
		
		initDic(Q, P, m);
	}
	
	lrs::~lrs() {
		lrs_free_dic(P, Q);
		lrs_free_dat(Q);
		lrs_close_quiet();
	}
	
	index_set lrs::allRatio(ind leave) {
		
		/* assign local variable to structures */
		matrix_t& A = P->A;
		ind* B = P->B;
		ind* Row = P->Row;
		ind* Col = P->Col;
		ind* minratio = Q->minratio;
		ind* inequality = Q->inequality;
		ind m = P->m;
		ind d = P->d;
		ind lastdv = Q->lastdv;
		
		ind cob, col;
		
		
		if ( (cob = findCob(leave)) < 0 ) 
			throw lrs_error("Failed to find cobasis for leaving index" + leave);
		
		col = Col[cob];
		
		ind degencount = 0;
		
		for (ind j = lastdv + 1; j <= m; j++) {
			/* search rows with negative coefficient in dictionary;
			 * minratio contains indices of min ratio cols */
			if ( negative(A[Row[j]][col]) ) minratio[degencount++] = j;
		}
		
		mpz_class Nmin, Dmin;
		
		ind ratiocol = 0;		/* column being checked, initially rhs */
		ind start = 0;			/* starting location in minratio array */
		ind bindex = d + 1;		/* index of next basic variable to consider */
		ind cindex = 0;			/* index of next cobasic variable to consider */
		ind basicindex = d;		/* index of basis inverse for current ratio 
								 * test, except d=rhs test */
		bool firstTime = true;	/* For ratio test, true on first pass, else 
								 * false */
		
		ind nstart = 0, ndegencount = 0;
		
		
		if ( B[bindex] == basicindex ) { 
			/* identity col in basis inverse */
			
			if ( minratio[start] == bindex ) {
				/* remove this index, all others stay */
				start++;
				degencount--;
			}
			bindex++;
			
		} else {
			/* perform ratio test on rhs or column of basis inverse */
			
			/* get next ratio column and increment cindex */
			for (ind j = start; j < start + degencount; j++) {
				ind i = Row[minratio[j]];	/* i is the row location of the 
											 * next basic variable */
				int comp = 1;	/* 1: lhs>rhs; 0: lhs=rhs; -1: lhs<rhs */
				
				if (firstTime) {
					firstTime = false; /* force new min ratio on first time */
				} else {
					if ( sgn(Nmin) > 0 || negative( A[i][ratiocol] ) ) {
						comp = ( sgn(Nmin) < 0 || positive( A[i][ratiocol] ) ) ?
							comprod( Nmin.get_mpz_t(), A[i][col], 
									 A[i][ratiocol], Dmin.get_mpz_t() )
							: -1;
					} else if ( sgn(Nmin) == 0 && zero( A[i][ratiocol] ) ) {
						comp = 0;
					}
					
					/* all signs reversed for rhs */
					if ( ratiocol == 0L ) comp = -comp;
				}
				
				if ( comp == 1 ) { /* new minimum ratio */
					nstart = j;
					Nmin = mpz_class( A[i][ratiocol] );
					Dmin = mpz_class( A[i][col] );
					ndegencount = 1;
				} else if ( comp == 0 ) { /* repeated minimum */
					minratio[nstart + ndegencount++] = minratio[j];
				}
			}
		}
		
		degencount = ndegencount;
		start = nstart;
		
		/* prepare return set */
		index_set rval(m+1);
		for (ind i = start; i < start + degencount; i++) {
			rval.set( inequality[ B[minratio[i]] - lastdv ] );
		}
		
		return rval;
	}
	
	ind lrs::findBas(ind enter) {
		ind j;
		ind lastdv = Q->lastdv;
		ind *inequality = Q->inequality;
		ind *B = P->B;
		ind m = P->m;
		
		for (j = lastdv+1; j <= m && inequality[ B[j]-lastdv ] != enter; j++);
		
		return (inequality[ B[j]-lastdv ] != enter) ? -1 : j;
	}
	
	ind lrs::findCob(ind leave) {
		ind j;
		ind lastdv = Q->lastdv;
		ind *inequality = Q->inequality;
		ind *C = P->C;
		ind d = P->d;
		
		for (j = 0; j < d && inequality[ C[j]-lastdv ] != leave; j++);
		
		return (inequality[ C[j]-lastdv ] != leave) ? -1 : j;
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
		
		index_set cobInd(m+1);
		index_set extraInc(m+1);
		
		for (i = 0; i < d; i++) {
			if (Col[i] == col) rflag = tempArray[i]; /* look for ray index */
			
			tempArray[i] = inequality[C[i] - lastdv];
			cobInd.set(tempArray[i]);
		}
		
		nIncidence = (col == 0) ? d : d-1;
		
		for (i = lastdv + 1; i <= m; i++) {
			if ( zero( A[Row[i]][0] ) ) {
				if ( (col == 0L) || zero( A[Row[i]][col] ) ) {
					extraInc.set(inequality[B[i] - lastdv]);
					nIncidence++;
				}
			}
		}
		
		mpz_class det(P->det);
		
		cobasis* cob = new cobasis(det, rflag, cobInd, nIncidence, extraInc);
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
		
		vector_mpz& output = *new vector_mpz(n);
		
		/* copy column 0 to output */
		copy( output[0], P->det );
		
		for (ind j = 1; j < n; j++) {
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
		
		return &output;
	}
	
	void lrs::initDat(lrs_dat* Q, ind n, ind d) {
		
		if (o.vRepresentation) {
			/* V-representation */
			Q->hull = true;
			Q->polytope = true;
			/* NOTE Symbal appears to always turn on geometric, LRS ignores it 
			 * for V-representation. I'm following Symbal here. */
			Q->geometric = true;
		} else {
			/* H-representation */
			Q->hull = false;
			
			/* copied from Symbal, equivalent to setting the geometric option 
			 * in LRS */
			Q->geometric = true;
		}
		
		Q->m = n; Q->n = d;
		
	}
	
	void lrs::initDic(lrs_dat* Q, lrs_dic* P, matrix const& mat) {
		
		ind m = Q->m;
		
		/* read matrix row by row */
		matrix& m_nc = const_cast<matrix &>(mat);
		for (ind i = 0; i < m; i++) {
			/* TODO allow for linearities */
			lrs_set_row_mp(
				P, Q, i+1, m_nc[i].num(), m_nc[i].den(), ge);
		}
		
	}
	
	ind lrs::lexRatio(ind leave) {
		ind cob;
		
		ind* Col = P->Col;
		ind* B = P->B;
		ind* inequality = Q->inequality;
		ind lastdv = Q->lastdv;
		
		if ( ( cob = findCob(leave) ) < 0 ) {
			return -1;
		}
		
		ind col = Col[cob];
		ind enter = ratio(P, Q, col);
		
		return ( enter > 0 ) ? inequality[ B[enter]-lastdv ] : -1;
	}
	
	void lrs::pivot(ind leave, ind enter) {
		
		ind cob = findCob(leave);
		if (cob < 0) throw lrs_error("Failed to find cobasis for pivot.");
		
		ind bas = findBas(enter);
		if (bas < 0) throw lrs_error("Failed to find basis for pivot.");
		
		::pivot(P, Q, bas, cob);
		::update(P, Q, &bas, &cob);
	}
	
	void lrs::printDict() {
		printA(P, Q);
	}
	
	void lrs::setCobasis(index_set& cob) {
		ind nlinearity = Q->nlinearity;
		ind* linearity = Q->linearity;
		ind* facet = Q->facet;
		ind m = Q->m;
		ind d = P->d;
		ind j = nlinearity;
		index_set_iter k = begin(cob);
		
		for ( ; j < d && k != end(cob); ++j, ++k) {
			facet[j] = *k;
			
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
