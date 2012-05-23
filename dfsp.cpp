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

#include "dfsp.hpp"
#include "fmt.hpp"
#include "gram.hpp"
#include "permUtils.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"

#include "permlib/permlib_api.h"
#include "permlib/permutation.h"

namespace basil {
	
	////////////////////////////////////////////////////////////////////
	//
	//  DFS Explorer Public Members
	//
	////////////////////////////////////////////////////////////////////
	
	dfsp::explorer::explorer(matrix& m, index_set& lin,
			permutation_group g, gram_matrix gram, dfsp_opts o)
			: l(m, lin, o.lrs_o), g(g), gramMat(gram), opts(o),
			  dim(m.dim()), rows(m.size()) {
		
		/* resize the cobasis cache to its proper size */
		cobasisCache.resize(opts.cacheSize);
		
		/* Default initialize remaining data members */
		basisOrbits = cobasis_map();
		cobasisGramMap = cobasis_gram_map();
		pathStack = pivot_stack();
		rayOrbits = coordinates_map();
		vertexOrbits = coordinates_map();
		vertexGramMap = vertex_gram_map();
		if ( opts.usesLocalStack ) {
			workStack = state_stack();
		}
	}
	
	bool dfsp::explorer::isKnownCobasis(
			cobasis_map& cobs, cobasis_gram_map& grams,
			index_set cob, vertex_data_ptr dat) {

		index_set_list possibleMatches =
				matchingCobasisInvariants(cobs, grams, cob, dat);
//int tid = omp_get_thread_num();
//#pragma omp critical(print)
//std::cout << "\t\t\t[" << tid << "] matchingInvariants: " << possibleMatches.size() << std::endl;

		/* if no known cobasis has invariants matching this one, it's new */
		if ( possibleMatches.size() == 0 ) return false;

		/* if a known cobasis (with matching invariants) is symmetric to this
		 * one, it's not new */
		if ( opts.assumesNoSymmetry || ! opts.stabSearch ) {
			/* no stabilizer search, so only one pass through the list */

			/* for each cobasis in the list to check for symmetry */
			for (index_set_list::iterator it = possibleMatches.begin();
					it != possibleMatches.end(); ++it) {

				/* the cobasis to check for symmetry */
				index_set old = *it;

				if ( cob == old ) {
					/* duplicate cobasis */
//#pragma omp critical(print)
//std::cout << "\t\t\t\t[" << tid << "] duplicate found" << std::endl;
					return true;
				}

				/* if no symmetry, check next cobasis */
				if ( opts.assumesNoSymmetry ) continue;

				/* look for a permutation in the permutation group that maps
				 * the incidence set of the cobasis we are trying to find to
				 * the incidence set of the known cobasis. */
				permutation_ptr act = permlib::setImage(g,
						plBegin(cob), plEnd(cob), plBegin(old), plEnd(old));

				/* This cobasis is symmetric to one we already know of */
				if (act) {
//#pragma omp critical(print)
//std::cout << "\t\t\t\t[" << tid << "] isomorph found" << std::endl;
					return true;
				}
			}
//#pragma omp critical(print)
//std::cout << "\t\t\t\t[" << tid << "] no isomorph found" << std::endl;
		} else {
			/* using stabilizer search, so need to check multiple grounds */

			/* represents the set [1..rows] */
			index_set allIndices = index_set(rows+1).set().set(0, false);

			/* for each possible size of superset of this cobasis */
			for (ind groundSize = cob.count() + 1; groundSize <= rows;
					groundSize++) {

				/* for each cobasis in the list to check for symmetry */
				for (index_set_list::iterator it = possibleMatches.begin();
						it != possibleMatches.end(); ++it) {

					/* the cobasis to check for symmetry */
					index_set old = *it;

					if ( cob == old ) {
						/* duplicate cobasis */
						return true;
					}

					/* Take the set union of the two cobases into ground */
					index_set ground = cob | old;
					/* Take the complement of ground into leftOut */
					index_set leftOut = allIndices - ground;

					while ( ground.count() < uind(groundSize) ) {
						/* take a random left out element and add to ground */
						ind randInd = lrs::pseudoRandomInd(leftOut);
						ground.set(randInd, true);
						leftOut.set(randInd, false);
					}

					/* get the set stabilizer of ground */
					permutation_group_ptr stab = permlib::setStabilizer(
							g, plBegin(ground), plEnd(ground));

					/* look for a permutation in the stabilizer group that maps
					 * the incidence set of the cobasis we are trying to find
					 * to the incidence set of the known cobasis. */
					permutation_ptr act = permlib::setImage(*stab,
							plBegin(cob), plEnd(cob),
							plBegin(old), plEnd(old));

					/* This cobasis is symmetric to one we already know of */
					if (act) return true;
				}
			}
		}

		/* if we can't find a matching cobasis, this one must be new */
		return false;
	}

	vertex_data_ptr dfsp::explorer::knownRay(
			coordinates_map& rays, vertex_data_ptr rep) {
		/* TODO think about including gram invariant here */
		
		/* incidence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (coordinates_map::iterator it = rays.begin(); 
				it != rays.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = it->second->inc;
			
			/* if we assume no symmetry, check for equal cobases */
			if ( opts.assumesNoSymmetry ) {
				if ( find == old ) return it->second; else continue;
			}
			
			/* Check incidence sets of equal size (a cheap invariant) */
			if ( find.count() != old.count() ) continue;
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the ray we are trying to find to the 
			 * incidence set of the known ray. */
			permutation_ptr act = permlib::setImage(
				g, plBegin(find), plEnd(find), 
				plBegin(old), plEnd(old));
			
			/* if such a permutation is found, return the known ray */
			if ( act ) return it->second;
		}
		
		/* no known ray that is equivalent up to symmetry */
		return vertex_data_ptr();
	}
	
	vertex_data_ptr dfsp::explorer::knownVertex(coordinates_map& verts,
			vertex_gram_map& grams, vertex_data_ptr rep) {
		coordinates_map::iterator found = verts.find(rep->coords);
		if ( found != verts.end() ) {
			/* duplicate vertex */
			return found->second;
		}

		/* if we assume no symmetry, it must be new */
		if ( opts.assumesNoSymmetry ) return vertex_data_ptr();

		/* List of vertices with matching invariants */
		vertex_data_list possibleMatches =
				matchingInvariants(verts, grams, rep);

		/* Test new by invariants */
		if ( possibleMatches.size() == 0 ) {
			return vertex_data_ptr();
		}

		/* Incidence set to find */
		index_set& find = rep->inc;

		/* for every known orbit representative that matches on invariants */
		for (vertex_data_list::const_iterator it = possibleMatches.begin();
				it != possibleMatches.end(); ++it) {

			/* incidence set to check */
			index_set& old = (*it)->inc;

			/* NOTE PermLib throws an assertion failure on setImage of
			 * differently sized sets, but the incidence set size is checked by
			 * matchingInvariants(). */

			/* look for a permutation in the global group that maps the
			 * incidence set of the vertex we are trying to find to the
			 * incidence set of the known vertex. */
			permutation_ptr act = permlib::setImage(
					g, plBegin(find), plEnd(find), plBegin(old), plEnd(old));
			/* if such a permuation is found, return the known vertex */
			if (act) {
				return *it;
			}
		}

		/* no known vertex that is equivalent up to symmetry */
		return vertex_data_ptr();
	}

	index_set_list dfsp::explorer::matchingCobasisInvariants(
			cobasis_map& cobs, cobasis_gram_map& grams,
			index_set cob, vertex_data_ptr dat) {

		/* list of cobases with matching invariants */
		index_set_list matches;

		if ( opts.gramVec ) {

			/* get set of cobases with matching gram vectors */
			cobasis_gram_range range =
					grams.equal_range( fastGramVec(gramMat, cob) );
//int tid = omp_get_thread_num();
//#pragma omp critical(print)
//std::cout << "\t\t\t[" << tid << "] has matching grams: " << ((range.second == range.first) ? "false" : "true") << std::endl;

			/* check invariants for each of these cobases */
			for (cobasis_gram_map::const_iterator it = range.first;
					it != range.second; ++it) {

				vertex_data_ptr oldDat = it->second.second;
//int matchInd = 0;

				if ( dat->inc.count() == oldDat->inc.count()
						&& dat->gram == oldDat->gram ) {
//#pragma omp critical(print)
//std::cout << "\t\t\t\t[" << tid << "] <" << matchInd << "> n = o = (" << dat->inc.count() << ", " << dat->gram << ")" << std::endl;
//matchInd++
					matches.push_back(it->second.first);
				} else {
//#pragma omp critical(print)
//std::cout << "\t\t\t\t[" << tid << "] <" << matchInd << "> n = (" << dat->inc.count() << ", " << dat->gram << "), o = (" << oldDat->inc.count() << ", " << oldDat->gram << ")" << std::endl;
//matchInd++
				}
			}

		} else {

			/* check invariants for each known cobasis */
			for (cobasis_map::const_iterator it = cobs.begin();
					it != cobs.end(); ++it) {

				vertex_data_ptr oldDat = it->second;

				if ( dat->inc.count() == oldDat->inc.count()
						&& dat->gram == oldDat->gram ) {
					matches.push_back(it->first);
				}

			}

		}

		return matches;
	}

	vertex_data_list dfsp::explorer::matchingInvariants(
			coordinates_map& verts, vertex_gram_map& grams,
			vertex_data_ptr rep) {

		/* list of vertices with matching invariants */
		vertex_data_list matches;

		if ( opts.gramVec ) {

			/* get set of cobases with matching gram vectors */
			vertex_gram_range range = grams.equal_range(rep->gram);

			/* check invariants for each of these cobases */
			for (vertex_gram_map::const_iterator it = range.first;
					it != range.second; ++it) {

				/* check invariants */
				if ( (*it->second).inc.count() == (*rep).inc.count() ) {
					matches.push_back(it->second);
				}
			}
		} else {

			/* check invariants for each known cobasis */
			for (coordinates_map::const_iterator it = verts.begin();
					it != verts.end(); ++it) {

				/* check invariants */
				if ( (*it->second).inc.count() == (*rep).inc.count() ) {
					matches.push_back(it->second);
				}
			}
		}

		return matches;
	}

	void dfsp::explorer::pivotTo(dfsp::pivot_stack const& target) {
		/* backtrack up state stack until it is only the prefix shared with the
		 * target stack, then fill the state stack back in with the remainder
		 * of the destination stack */

		//first make sure the size is the same
		while ( pathStack.size() > target.size() ) {
			pivot const& p = pathStack.back();
			l.pivot(p.enter, p.leave);
			pathStack.pop_back();
		}

		//then go back to the shared prefix
		ind i = pathStack.size() - 1;
		while ( i >= 0 ) {
			pivot const& p = pathStack.back();
			l.pivot(p.enter, p.leave);

			//check for last cobasis before shared prefix
			if ( p.cob == target.at(i).cob ) {
				pathStack.pop_back();
				--i;
				break;
			}

			pathStack.pop_back();
			--i;
		}

		//then fill the state stack back in
		while ( ++i < (ind)target.size() ) {
			pivot const& p = target.at(i);
			pathStack.push_back(p);
			l.pivot(p.leave, p.enter);
		}
	}

	////////////////////////////////////////////////////////////////////
	//
	//  DFS Public Members
	//
	////////////////////////////////////////////////////////////////////
	
	dfsp::dfsp(matrix& m, index_set& lin, permutation_group& g,
			gram_matrix& gram, dfsp_opts o) : globalM(m), globalLin(lin),
			globalG(g), globalOpts(o),
			globalDim(m.dim()), globalRows(m.size()), globalGramMat(gram) {
		
		/* set up algorithm globals */
		initGlobals();
	}
	
	bool dfsp::doDfs() {
		
		/* set algorithm start time */
		start_time = std::clock();
#ifdef BAS_WALLTIME
		gettimeofday(&wall_start_time, 0);
#endif /* BAS_WALLTIME */
		
		/* Algorithm success */
		int nThreads = 0;
		int nSuccess = 0;
		/* Number of stalled threads */
		int nWaiting = 0;
		
		#pragma omp parallel reduction(+:nSuccess)
		{
		#pragma omp master
		{
		nThreads = omp_get_num_threads();
		} /* omp master */
		
		/* Set up thread locals */
		explorer ex(globalM, globalLin, globalG, globalGramMat, globalOpts);
		
		/* print initial dictionary */
		#pragma omp master
		{
		if ( globalOpts.showsAllDicts ) ex.l.printDict();
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
		if ( globalOpts.showsAllDicts ) ex.l.printDict();
		} /* omp master */
		
		/* go to the initial cobasis, if supplied */
		if ( globalOpts.firstCobasis ) {
			ex.l.setCobasis( *globalOpts.firstCobasis );
		}
		
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
		
		initialCobasis = cob->cob;
		ex.cobasisCache.insert(initialCobasis);

		#pragma omp master
		{
		if ( globalOpts.printTrace ) {
			globalOpts.output() << "#I initial basis: " << fmt( cob->cob )
					<< " " << *sol << "\n";
		}

		addVertex(dat);
		getRays(ex);
		} /* omp master */
		
		/* DFS the edge graph */
		pushNewEdges(ex, cob->cob);

		/* construct empty pivot */
		pivot_stack ps;
		bool waiting = false;
		bool working = true;

		while ( working ) {
			if ( globalOpts.usesLocalStack && ! ex.workStack.empty() ) {
				/* pop the pivot to explore off of the local stack */
				ps = ex.workStack.back(); ex.workStack.pop_back();
				/* share local work stack */
				if ( ! ex.workStack.empty() ) {
					#pragma omp critical(stacks)
					{
						if ( nWaiting > 0 ) {
							/* unload local work stack into global if needed */
							do {
								globalWorkStack.push_back(ex.workStack.back());
								ex.workStack.pop_back();
							} while ( ! ex.workStack.empty() );
						}
					} /* omp critical(stacks) */
				}
				/* decide whether to run again or not */
				working = ( globalOpts.basisLimit > ex.basisOrbits.size() );
			} else {
				#pragma omp critical(stacks)
				{
				/* check for work */
				if ( globalWorkStack.empty() ) {
					if ( ! waiting ) {
						waiting = true;
						++nWaiting;
					}
				} else {
					/* pop the pivot to the edge to explore off the global
					 * stack */
					ps = globalWorkStack.back(); globalWorkStack.pop_back();
					if ( waiting ) {
						waiting = false;
						--nWaiting;
					}
				}
				/* decide whether to run again or not */
				working = ( nWaiting < nThreads
						&& globalOpts.basisLimit > ex.basisOrbits.size() );
				} /* omp critical(stacks) */
			}

			/* poll for new work if we don't have any */
			if ( waiting ) continue;

			/* pivot explorer to the new basis */
			ex.pivotTo(ps);

			/* get the last pivot and current cobasis */
			pivot& p = ex.pathStack.back();

			/* print the dictionary pivoted to */
			if ( globalOpts.showsAllDicts ) {
				#pragma omp critical(print)
				{
				ex.l.printDict();
				} /* omp critical(print) */
			}
			if ( globalOpts.printTrace ) {
				#pragma omp critical(print)
				{
				globalOpts.output() << "#I traversing " << fmt( p.cob )
							<< " through (" << p.leave << "," << p.enter
							<< ")\n";
				} /* omp critical(print) */
			}

			/* get the new cobasis */
			cobasis_ptr cob(ex.l.getCobasis(0));
			/* get the new rays */
			getRays(ex);

			/* Add new vertex representations / cobases adjacent to the new
			 * vertex to the work stack */
			pushNewEdges(ex, cob->cob);
		}

		nSuccess += ( globalOpts.basisLimit >= ex.basisOrbits.size() );
		} /* omp parallel */
		
		/* set algorithm end time */
		diff_time = std::clock() - start_time;
#ifdef BAS_WALLTIME
		gettimeofday(&wall_end_time, 0);
#endif /* BAS_WALLTIME */
		
		return ( nThreads == nSuccess );
	}
	
	////////////////////////////////////////////////////////////////////
	// Query methods for after completion of doDfs()
	////////////////////////////////////////////////////////////////////
	
	cobasis_map dfsp::getBasisOrbits() const {
		return cobasis_map(globalBasisOrbits.begin(), 
				globalBasisOrbits.end());
	}
	
	ind dfsp::getDimension() const { return globalDim - 1; }
	
	index_set dfsp::getInitialCobasis() const { return initialCobasis; }
	
	bool dfsp::isFinished() const { return !hitMaxBasis; }
	
	coordinates_map dfsp::getRayOrbits() const {
		return coordinates_map(globalRayOrbits.begin(), 
				globalRayOrbits.end());
	}
	
	std::clock_t dfsp::getRunningTime() const
		{ return diff_time / clocks_per_ms; }
	
#ifdef BAS_WALLTIME
	long dfsp::getWallTime() const {
		double d_start =
				wall_start_time.tv_sec * 1000000 + wall_start_time.tv_usec;
		double d_end =
				wall_end_time.tv_sec * 1000000 + wall_end_time.tv_usec;
		return (long)((d_end - d_start)/1000);
	}
#endif /* BAS_WALLTIME */

	permutation_group const& dfsp::getSymmetryGroup() const { return globalG; }
	
	coordinates_map dfsp::getVertexOrbits() const {
		return coordinates_map(globalVertexOrbits.begin(), 
				globalVertexOrbits.end());
	}
	
	gram_matrix const& dfsp::getGramMat() const { return globalGramMat; }
	
	
	////////////////////////////////////////////////////////////////////
	//
	//  DFS Private Members
	//
	////////////////////////////////////////////////////////////////////
	
	void dfsp::addCobasis(index_set const& cob, vertex_data_ptr dat) {
		
		unsigned int oSize;
		#pragma omp critical(globals)
		{
		globalBasisOrbits.push_back(std::make_pair(cob, dat));
		oSize = globalBasisOrbits.size();
		} /* omp critical(globals) */
		
		/* print cobasis, if option set */
		if ( globalOpts.printBasis && oSize % globalOpts.printBasis == 0 ) {
			#pragma omp critical(print)
			{
			std::ostream& out = globalOpts.output();
			out << "# cobases: " << oSize << " (" << currentTime() 
					<< " ms)";
			if ( globalOpts.printNew ) {
				out << " " << fmt( cob );
				if ( globalOpts.debugGram ) out << " " << dat->gram;
			}
			out << std::endl;
			} /* omp critical(print) */
		}
	}
	
	void dfsp::addVertex(vertex_data_ptr dat) {
		
		unsigned int oSize;
		#pragma omp critical(globals)
		{
		/* map the rationalization of the coordinates to the vertex 
		 * data */
		globalVertexOrbits.push_back(std::make_pair(dat->coords, dat));
		oSize = globalVertexOrbits.size();
		} /* omp critical(globals) */
		
		/* for each defined cobasis, map it to the vertex data */
		for (std::set<index_set>::iterator it = dat->cobs.begin();
				it != dat->cobs.end(); ++it) {
			addCobasis(*it, dat);
		}
		
		/* print vertex, if option set */
		if ( globalOpts.printVertex && oSize % globalOpts.printVertex == 0 ) {
			#pragma omp critical(print)
			{
			std::ostream& out = globalOpts.output();
			out << "# vertices: " << oSize << " (" << currentTime() 
					<< " ms)";
			if ( globalOpts.printNew ) {
				out << " " << dat->coords;
				if ( globalOpts.debugGram ) out << " " << dat->gram;
			}
			out << std::endl;
			} /* omp critical(print) */
		}
	}
	
	void dfsp::getRays(dfsp::explorer& ex) {
		for (ind j = 1; j <= realDim; j++) {
			vector_mpz_ptr s( ex.l.getSolution(j) );
			
			if (s) {
				cobasis_ptr c( ex.l.getCobasis(j) );
				vertex_data_ptr dat( rayData(c, s) );
				
				/* Not a new ray, by local cache */
				if ( ex.knownRay(ex.rayOrbits, dat) ) continue;
				
				while ( true ) {
					coordinates_map newRayOrbits;
					uind oSize;

					#pragma omp critical(rays)
					{
					if ( ex.rayOrbits.size() == globalRayOrbits.size() ) {
						/* local up to date, add to global */
						globalRayOrbits.push_back(
								std::make_pair(dat->coords, dat));
					} else {
						/* update local to current global state */
						newRayOrbits.insert(
								globalRayOrbits.begin() + ex.rayOrbits.size(),
								globalRayOrbits.end());
					}
					oSize = globalRayOrbits.size();
					} /* omp critical rays */

					if ( newRayOrbits.empty() ) {
						/* found a new ray */

						ex.rayOrbits.insert(std::make_pair(dat->coords,
								boost::make_shared<vertex_data>(dat->coords,
										dat->inc, dat->cobs, dat->det,
										dat->gram)));

						if ( globalOpts.printRay
								&& oSize % globalOpts.printRay == 0 ) {
							#pragma omp critical(print)
							{
							std::ostream& out = globalOpts.output();
							out << "# rays: " << oSize << " ("
									<< currentTime() << " ms)";
							if ( globalOpts.printNew ) {
								out << " " << dat->coords;
								if ( globalOpts.debugGram ) {
									out << " " << dat->gram;
								}
							}
							out << std::endl;
							} /* omp critical print */
						}

						break;
					} else {
						/* update locals with new globals */

						for (coordinates_map::iterator it =
								newRayOrbits.begin(); it != newRayOrbits.end();
								++it) {
							vertex_data& val = *(it->second);

							ex.rayOrbits.insert(std::make_pair(val.coords,
									boost::make_shared<vertex_data>(val.coords,
											val.inc, val.cobs, val.det,
											val.gram)));
						}

						/* Not a new ray, by new local cache */
						if ( ex.knownRay(newRayOrbits, dat) ) break;
					}
				}
			}
		}
	}
	
	void dfsp::initGlobals() {
		/* account for flipability of arrangement gram matrix */
		if ( globalOpts.aRepresentation ) globalGramMat = globalGramMat.abs();
		
		/* Default initialize remaining data members */
		globalBasisOrbits = cobasis_list();
		hitMaxBasis = false;
		initialCobasis = index_set();
		globalRayOrbits = coordinates_list();
		realDim = 0;
		globalVertexOrbits = coordinates_list();
		globalWorkStack = state_stack();
	}
	
	bool dfsp::knownOrAddNewCobasis(dfsp::explorer& ex,
			index_set cob, vertex_data_ptr dat) {

//int tid = omp_get_thread_num();
//#pragma omp critical(print)
//std::cout << "\t[" << tid << "] exploring cobasis " << fmt( cob ) << std::endl;
		/* check cobasis against local store */
		bool known =
				ex.isKnownCobasis(ex.basisOrbits, ex.cobasisGramMap, cob, dat);
//if ( known ) {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] rejected by locals: l = " << ex.basisOrbits.size() << std::endl;
//}

		while ( ! known ) {
			cobasis_map newBasisOrbits;
			cobasis_gram_map newCobasisGramMap;

			#pragma omp critical(cobases)
			{
			if ( ex.basisOrbits.size() == globalBasisOrbits.size() ) {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] accepted: l = g = " << ex.basisOrbits.size() << std::endl;
				/* local up to date, add to global */
				addCobasis(cob, dat);
			} else {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] delayed: l = " << ex.basisOrbits.size() << ", g = " << globalBasisOrbits.size() << std::endl;
				/* update local to current global state */
				newBasisOrbits.insert(
						globalBasisOrbits.begin() + ex.basisOrbits.size(),
						globalBasisOrbits.end());
			}
			} /* omp critical(cobases) */

			if ( newBasisOrbits.empty() ) {
				/* found a new cobasis */

				vertex_data_ptr newDat = boost::make_shared<vertex_data>(
						dat->coords, dat->inc, dat->cobs, dat->det, dat->gram);

				ex.basisOrbits.insert(std::make_pair(cob, newDat));
				if ( globalOpts.gramVec ) {
//#pragma omp critical(print)
//std::cout << "\t\t\t[" << tid << "] inserted " << fmt( cob ) <<  fastGramVec(globalGramMat, cob) << std::endl;
					ex.cobasisGramMap.insert(std::make_pair(
							fastGramVec(globalGramMat, cob),
							std::make_pair(cob, newDat)));
				}

				return true;
			} else {
				/* update locals with new globals */

				for (cobasis_map::iterator it = newBasisOrbits.begin();
						it != newBasisOrbits.end(); ++it) {
					const index_set& key = it->first;
					vertex_data& val = *(it->second);
					vertex_data_ptr newDat = boost::make_shared<vertex_data>(
							val.coords, val.inc, val.cobs, val.det, val.gram);

					ex.basisOrbits.insert(std::make_pair(key, newDat));
					if ( globalOpts.gramVec ) {
						gram_matrix gram = fastGramVec(globalGramMat, key);
//#pragma omp critical(print)
//std::cout << "\t\t\t[" << tid << "] updated " << fmt( key ) <<  gram << std::endl;
						ex.cobasisGramMap.insert(
								std::make_pair(gram,
										std::make_pair(key, newDat)));
						newCobasisGramMap.insert(
								std::make_pair(gram,
										std::make_pair(key, newDat)));
					}
				}

				/* test vertex against new local cache */
				known = ex.isKnownCobasis(
						newBasisOrbits, newCobasisGramMap, cob, dat);
//if ( known ) {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] rejected by updates: u = " << newBasisOrbits.size() << ", l = " << ex.basisOrbits.size() << std::endl;
//}
			}
		}

		return false;
	}

	dfsp::vertex_data_known dfsp::knownOrAddNewVertex(dfsp::explorer& ex,
			vertex_data_ptr rep) {
//int tid = omp_get_thread_num();
//#pragma omp critical(print)
//std::cout << "\t[" << tid << "] exploring vertex " << rep->coords << std::endl;
		/* check vertex against local store */
		vertex_data_ptr known =
				ex.knownVertex(ex.vertexOrbits, ex.vertexGramMap, rep);
//if ( known ) {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] rejected by locals: l = " << ex.vertexOrbits.size() << std::endl;
//}

		while ( ! known ) {
			coordinates_map newVertexOrbits;
			vertex_gram_map newVertexGramMap;

			#pragma omp critical(vertices)
			{
			if ( ex.vertexOrbits.size() == globalVertexOrbits.size() ) {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] accepted: l = g = " << ex.vertexOrbits.size() << std::endl;
				/* local up to date, add to global */
				addVertex(rep);
			} else {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] delayed: l = " << ex.vertexOrbits.size() << ", g = " << globalVertexOrbits.size() << std::endl;
				/* update local to current global state */
				newVertexOrbits.insert(
						globalVertexOrbits.begin() + ex.vertexOrbits.size(),
						globalVertexOrbits.end());
			}
			} /* omp critical(vertices) */

			if ( newVertexOrbits.empty() ) {
				/* found a new vertex */

				vertex_data_ptr dat = boost::make_shared<vertex_data>(
						rep->coords, rep->inc, rep->cobs, rep->det, rep->gram);
				ex.vertexOrbits.insert(std::make_pair(dat->coords, dat));
				if ( globalOpts.gramVec ) {
					ex.vertexGramMap.insert(std::make_pair(dat->gram, dat));
				}

				return vertex_data_known(dat, true);
			} else {
				/* update locals with new globals */

				for (coordinates_map::iterator it = newVertexOrbits.begin();
						it != newVertexOrbits.end(); ++it) {
					vertex_data& val = *(it->second);
					vertex_data_ptr dat = boost::make_shared<vertex_data>(
							val.coords, val.inc, val.cobs, val.det, val.gram);

					ex.vertexOrbits.insert(std::make_pair(dat->coords, dat));
					if ( globalOpts.gramVec ) {
						ex.vertexGramMap.insert(std::make_pair(dat->gram, dat));
						newVertexGramMap.insert(std::make_pair(dat->gram, dat));
					}
				}

				/* test vertex against new local cache */
				known = ex.knownVertex(newVertexOrbits, newVertexGramMap, rep);
//if ( known ) {
//#pragma omp critical(print)
//std::cout << "\t\t[" << tid << "] rejected by updates: u = " << newVertexOrbits.size() << ", l = " << ex.vertexOrbits.size() << std::endl;
//}
			}
		}

		return vertex_data_known(known, false);
	}

	void dfsp::pushNewEdges(explorer& ex, index_set& oldCob) {

		/* for each index in the old cobasis */
		for (index_set_iter it = lrs::begin(oldCob);
				it != lrs::end(oldCob); ++it) {

			/* the leaving index */
			ind leave = *it;
			/* the appropriate entering indices */
			index_set entering(oldCob.size());
			/* the entering index */
			ind enter;

			if ( globalOpts.aRepresentation ) {
				/* use arrangement pivot selection */
				entering = ex.l.arrangementRatio(leave);
			} else if ( globalOpts.lexOnly ) {
				/* calculate entering index lexicographically (BAD) */
				enter = ex.l.lexRatio(leave);
				if (enter >= 0) entering.set(enter); else continue;
			} else {
				/* calculate set of valid entering indices */
				entering = ex.l.allRatio(leave);
			}

			if ( globalOpts.printTrace ) {
				#pragma omp critical(print)
				{
				globalOpts.output() << "#I for leaving index { " << leave
							<< " } possible entering " << fmt( entering )
							<< "\n";
				} /* omp critical(print) */
			}

			/* for each valid entering index */
			for (index_set_iter jt = lrs::begin(entering);
					jt != lrs::end(entering); ++jt) {

				enter = *jt;

				/* Do the given pivot, then get the cobasis for the new edge */
				ex.l.pivot(leave, enter);
				cobasis_ptr cob(ex.l.getCobasis(0));
				vector_mpz_ptr sol(ex.l.getVertex());
				if ( globalOpts.showsAllDicts ) {
					#pragma omp critical(print)
					{
					globalOpts.output() << "\nPivot: " << leave << "=>"
							<< enter;
					ex.l.printDict();
					} /* omp critical(print) */
				}
				/* pivot back */
				ex.l.pivot(enter, leave);

				/* avoid expensive symmetry calculations by caching recently
				 * seen cobases.
				 *
				 * As the time to insert or lookup should be similar, using
				 * insert here instead of lookup saves a call.
				 */
				if ( ! ex.cobasisCache.insert(cob->cob) ) {

					/* if this cobasis is not in the cache, add it */
					/* ex.cobasisCache.insert(cob->cob); */

					/* calculate invariants of new cobasis */
					vertex_data_ptr dat(vertexData(cob, sol));
					vertex_data_known vert = knownOrAddNewVertex(ex, dat);
					if ( vert.second ) {
						/* a new vertex */

						pivot_stack newWork = ex.pathStack;
						newWork.push_back(pivot(oldCob, leave, enter));
						if ( globalOpts.usesLocalStack ) {
							ex.workStack.push_back(newWork);
						} else {
							#pragma omp critical(stacks)
							{
							globalWorkStack.push_back(newWork);
							} /* omp critical(stacks) */
						}

						if ( globalOpts.printTrace ) {
							#pragma omp critical(print)
							{
							globalOpts.output() << "#I pushing new vertex: "
									<< fmt( cob->cob ) << " " << *sol << "\n";
							} /* omp critical(print) */
						}
					} else if ( dat->coords == vert.first->coords
							|| ! globalOpts.dualFacetTrick ) {

						/* if this is a new cobasis for a previously seen
						 * vertex, and we are not employing the dual facet
						 * trick to prune the search tree, add the cobasis to
						 * the search stack if it is unique */

						if ( knownOrAddNewCobasis(ex, cob->cob, vert.first) ) {

							pivot_stack newWork = ex.pathStack;
							newWork.push_back(pivot(oldCob, leave, enter));
							if ( globalOpts.usesLocalStack ) {
								ex.workStack.push_back(newWork);
							} else {
								#pragma omp critical(stacks)
								{
								globalWorkStack.push_back(newWork);
								} /* omp critical(stacks) */
							}

							if ( globalOpts.printTrace ) {
								#pragma omp critical(print)
								{
								globalOpts.output()
										<< "#I pushing new cobasis: "
										<< fmt( cob->cob ) << " " << *sol
										<< "\n";
								} /* omp critical(print) */
							}
						}
					} else {

						/* we assume that if the new cobasis is defining a
						 * different vertex, but that vertex is symmetric, then
						 * its neighbours will be symmetric to those of the
						 * known vertex. Prune via the dual facet trick. */
						if ( globalOpts.printTrace ) {
							#pragma omp critical(print)
							{
							globalOpts.output() << "#I ignoring cobasis "
									<< fmt( cob->cob ) << " by dual facet "
									"trick\n";
							} /* omp critical(print) */
						}
					}
				} else if ( globalOpts.printTrace ) {
					#pragma omp critical(print)
					{
					globalOpts.output() << "#I seen cobasis "
							<< fmt( cob->cob ) << " before\n";
					} /* omp critical(print) */
				}
			}
		}
	}

	vertex_data_ptr dfsp::rayData(dfsp::cobasis_ptr cob,
			dfsp::vector_mpz_ptr coords) {
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
	
	vertex_data_ptr dfsp::vertexData(dfsp::cobasis_ptr cob,
			dfsp::vector_mpz_ptr coords) {
		
		/* union of the cobasis and extra incidence of the cobasis 
		 * data */
		index_set inc = cob->cob | cob->extraInc;
		/* gram matrix, or empty if option off */
		gram_matrix gram = ( globalOpts.gramVec ) ?
				fastGramVec(globalGramMat, inc) : gram_matrix();
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				coords->rationalization(), inc, cob->cob, 
				abs(cob->det), gram);
		
		return dat;
	}

	gram_matrix fastGramVec(gram_matrix& gramMat, index_set inc) {
		/* restrict the inner product matrix to the incidence set, then
		 * sort it to the canonical representation of that matrix */
		return gramMat.restriction(inc).sort();
	}

} /* namespace basil */
