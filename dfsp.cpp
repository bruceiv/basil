#include <algorithm>
#include <ctime>
#include <deque>
#include <iterator>
#include <ostream>
#include <set>
#include <vector>

#include <boost/make_shared.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include <omp.h>

#include <permlib/permlib_api.h>
#include <permlib/permutation.h>

#include "dfsp.hpp"
#include "fmt.hpp"
#include "gram.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"

namespace basil {
	
	////////////////////////////////////////////////////////////////////
	//
	//  DFS Explorer Public Members
	//
	////////////////////////////////////////////////////////////////////
	
	dfs::explorer::explorer(matrix& m, index_set& lin, 
			permutation_group& g, gram_matrix& gram, dfs_opts o) 
			: l(m, lin, o.lrs_o), g(g), gramMat(gram) {
		
	}
	
	////////////////////////////////////////////////////////////////////
	//
	//  DFS Public Members
	//
	////////////////////////////////////////////////////////////////////
	
	dfs::dfs(matrix& m, index_set& lin, permutation_group& g, 
			gram_matrix& gram, dfs_opts o) : m(m), lin(lin), g(g), 
			opts(o), dim(m.dim()), rows(m.size()), gramMat(gram) { 
		
		/* set up algorithm globals */
		initGlobals();
	}
	
	bool dfs::doDfs() {
		
		/* set algorithm start time */
		start_time = std::clock();
		
		/* Algorithm success */
		int nThreads = 0;
		int nSuccess = 0;
		
		#pragma omp parallel reduction(+:nSuccess)
		{
		#pragma omp master
		{
		nThreads = omp_get_num_threads();
		} /* omp master */
		
		/* Set up thread locals */
		explorer ex(m, lin, g, gramMat, opts);
		
		/* print initial dictionary */
		#pragma omp master
		{
		if (opts.showsAllDicts) ex.l.printDict();
		} /* omp master */ 
		
		////////////////////////////////////////////////////////////////
		// get initial cobasis / vertex either from options or from 
		// LRS.
		// 
		// NOTE: this is a bit kludgy, due to the difficulty of copying 
		// LRS tableaux, so each thread duplicates the work of finding 
		// an initial basis done by the master.
		////////////////////////////////////////////////////////////////
		
		/* get initial solution from LRS */
		if ( ! ex.l.getFirstBasis() ) 
			throw dfs_error("LRS failed to find first basis.");
		
		/* print initial solution */
		#pragma omp master
		{
		if (opts.showsAllDicts) ex.l.printDict();
		} /* omp master */
		
		/* go to the initial cobasis, if supplied */
		if ( opts.firstCobasis ) ex.l.setCobasis( *opts.firstCobasis );
		
		/* get true problem dimension */
		#pragma omp master
		{
		realDim = ex.l.getRealDim();
		} /* omp master */
		
		/* add this vertex / cobasis as a representative of its orbit 
		 * (obviously there aren't any others yet) */
		cobasis_ptr cob(ex.l.getCobasis(0));
		vector_mpz_ptr sol(ex.l.getVertex());
		vertex_data_ptr dat(vertexData(cob, sol));
		
		#pragma omp master
		{
		if ( opts.printTrace ) {
			opts.output() << "#I initial basis: " << fmt( cob->cob ) 
					<< " " << *sol << "\n";
		}
		
		initialCobasis = cob->cob;
		
		addVertex(dat);
		
		//TODO FIXME getRays(ex);
		
		cobasisCache.insert(initialCobasis);
		} /* omp master */
		
		/* DFS the edge graph, returning whether it successfully completes */
		//TODO FIXME res = dfsFromRoot(cob);
		nSuccess += true;
		} /* omp parallel */
		
		/* set algorithm end time */
		diff_time = std::clock() - start_time;
		
		return ( nThreads == nSuccess );
	}
	
	////////////////////////////////////////////////////////////////////
	// Query methods for after completion of doDfs()
	////////////////////////////////////////////////////////////////////
	
	dfs::cobasis_map const& dfs::getBasisOrbits() const 
		{ return basisOrbits; }
	
	ind dfs::getDimension() const { return dim - 1; }
	
	index_set dfs::getInitialCobasis() const { return initialCobasis; }
	
	bool dfs::isFinished() const { return !hitMaxBasis; }
	
	dfs::coordinates_map const& dfs::getRayOrbits() const 
		{ return rayOrbits; }
	
	std::clock_t dfs::getRunningTime() const 
		{ return diff_time / clocks_per_ms; }
	
	permutation_group const& dfs::getSymmetryGroup() const { return g; }
	
	dfs::coordinates_map const& dfs::getVertexOrbits() const 
		{ return vertexOrbits; }
	
	gram_matrix const& dfs::getGramMat() const { return gramMat; }
	
	dfs::cobasis_gram_map const& dfs::getCobasisGramMap() const 
		{ return cobasisGramMap; }
	
	dfs::vertex_gram_map const& dfs::getVertexGramMap() const 
		{ return vertexGramMap; }
	
	
	////////////////////////////////////////////////////////////////////
	//
	//  DFS Private Members
	//
	////////////////////////////////////////////////////////////////////
	
	void dfs::addCobasis(index_set const& cob, 
			dfs::vertex_data_ptr dat) {
		
		unsigned int oSize;
		#pragma omp critical(globals)
		{
		basisOrbits.insert(std::make_pair(cob, dat));
		oSize = basisOrbits.size();
		
		/* add gram vector, if option set */
		if (opts.gramVec) {
			cobasisGramMap.insert(
				std::make_pair(
					fastGramVec(cob), std::make_pair(cob, dat))
			);
		}
		} /* omp critical(globals) */
		
		/* print cobasis, if option set */
		if ( opts.printBasis && oSize % opts.printBasis == 0 ) {
			#pragma omp critical(print)
			{
			std::ostream& out = opts.output();
			out << "# cobases: " << oSize << " (" << currentTime() 
					<< " ms)";
			if ( opts.printNew ) {
				out << " " << fmt( cob );
				if ( opts.debugGram ) out << " " << dat->gram;
			}
			out << std::endl;
			} /* omp critical(print) */
		}
	}
	
	void dfs::addVertex(dfs::vertex_data_ptr dat) {
		
		unsigned int oSize;
		#pragma omp critical(globals)
		{
		/* map the rationalization of the coordinates to the vertex 
		 * data */
		vertexOrbits.insert(std::make_pair(dat->coords, dat));
		oSize = vertexOrbits.size();
		
		/* add gram vector, if option set */
		if (opts.gramVec) {
			vertexGramMap.insert(
					std::make_pair(fastGramVec(dat->inc), dat));
		}
		} /* omp critical(globals) */
		
		/* for each defined cobasis, map it to the vertex data */
		for (std::set<index_set>::iterator it = dat->cobs.begin();
				it != dat->cobs.end(); ++it) {
			addCobasis(*it, dat);
		}
		
		/* print vertex, if option set */
		if ( opts.printVertex && oSize % opts.printVertex == 0 ) {
			#pragma omp critical(print)
			{
			std::ostream& out = opts.output();
			out << "# vertices: " << oSize << " (" << currentTime() 
					<< " ms)";
			if ( opts.printNew ) { 
				out << " " << dat->coords;
				if ( opts.debugGram ) out << " " << dat->gram;
			}
			out << std::endl;
			} /* omp critical(print) */
		}
	}
	
	gram_matrix dfs::fastGramVec(index_set inc) {
		/* restrict the inner product matrix to the incidence set, then 
		 * sort it to the canonical representation of that matrix */
		return gramMat.restriction(inc).sort();
	}
	
	/* TODO rework me
	void dfs::getRays(dfs::explorer ex) {
		for (ind j = 1; j <= realDim; j++) {
			vector_mpz_ptr s( ex.l.getSolution(j) );
			
			if (s) {
				cobasis_ptr c( ex.l.getCobasis(j) );
				vertex_data_ptr dat( rayData(c, s) );
				if ( ! knownRay(dat) ) {
					rayOrbits.insert(std::make_pair(dat->coords, dat));
					
					if ( opts.printRay 
							&& rayOrbits.size() % opts.printRay == 0 ) {
						std::ostream& out = opts.output();
						out << "# rays: " << rayOrbits.size() << " (" 
								<< currentTime() << " ms)";
						if ( opts.printNew ) {
							out << " " << dat->coords;
							if ( opts.debugGram ) out << " " << dat->gram;
						}
						out << std::endl;
					}
				}
			}
		}
	}
	*/
	
	void dfs::initGlobals() {
		/* resize the cobasis cache to its proper size */
		cobasisCache.resize(opts.cacheSize);
		/* account for flipability of arrangement gram matrix */
		if ( opts.aRepresentation ) gramMat = gramMat.abs();
		
		/* Default initialize remaining data members */
		basisOrbits = cobasis_map();
		cobasisGramMap = cobasis_gram_map();
		hitMaxBasis = false;
		initialCobasis = index_set();
		rayOrbits = coordinates_map();
		realDim = 0;
		vertexOrbits = coordinates_map();
		vertexGramMap = vertex_gram_map();
	}
	
	dfs::vertex_data_ptr dfs::rayData(dfs::cobasis_ptr cob, 
			dfs::vector_mpz_ptr coords) {
		/* TODO look at including gramVec from dfs.gap RayRep() */
		
		/* union of the cobasis and extra incidence of the cobasis 
		 * data */
		index_set inc = cob->cob | cob->extraInc;
		/* less the ray index */
		inc.set(cob->ray, false);
		/* ignore gram vector for the moment */
		gram_matrix gram = gram_matrix();
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				coordinates(*coords), inc, cob->cob, abs(cob->det), 
				gram);
		
		return dat;
	}
	
	dfs::vertex_data_ptr dfs::vertexData(dfs::cobasis_ptr cob, 
			dfs::vector_mpz_ptr coords) {
		
		/* union of the cobasis and extra incidence of the cobasis 
		 * data */
		index_set inc = cob->cob | cob->extraInc;
		/* gram matrix, or empty if option off */
		gram_matrix gram = ( opts.gramVec ) ? 
				fastGramVec(inc) : gram_matrix();
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				coords->rationalization(), inc, cob->cob, 
				abs(cob->det), gram);
		
		return dat;
	}

} /* namespace basil */
