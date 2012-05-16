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
#include "permUtils.hpp"

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
			permutation_group g, gram_matrix gram, dfs_opts o)
			: l(m, lin, o.lrs_o), g(g), gramMat(gram), opts(o),
			  dim(m.dim()), rows(m.size()) {
		
		/* resize the cobasis cache to its proper size */
		cobasisCache.resize(opts.cacheSize);
		
		/* Default initialize remaining data members */
		basisOrbits = cobasis_map();
		cobasisGramMap = cobasis_gram_map();
		cobasisUpdate = 0;
		rayOrbits = coordinates_map();
		rayUpdate = 0;
		vertexOrbits = coordinates_map();
		vertexGramMap = vertex_gram_map();
		vertexUpdate = 0;
	}
	
	bool dfs::explorer::isKnownCobasis(
			dfs::cobasis_map cobs, dfs::cobasis_gram_map grams,
			index_set cob, dfs::vertex_data_ptr dat) {

		index_set_list possibleMatches =
				matchingCobasisInvariants(cobs, grams, cob, dat);

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
				if (act) return true;
			}
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

	dfs::vertex_data_ptr dfs::explorer::knownRay(
			dfs::coordinates_map rays, dfs::vertex_data_ptr rep) {
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
	
	dfs::vertex_data_ptr dfs::explorer::knownVertex(dfs::coordinates_map verts,
			dfs::vertex_gram_map grams, dfs::vertex_data_ptr rep) {

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
		if ( possibleMatches.size() == 0 ) return vertex_data_ptr();

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
			if (act) return *it;
		}

		/* no known vertex that is equivalent up to symmetry */
		return vertex_data_ptr();
	}

	dfs::index_set_list dfs::explorer::matchingCobasisInvariants(
			dfs::cobasis_map cobs, dfs::cobasis_gram_map grams,
			index_set cob, dfs::vertex_data_ptr dat) {

		/* list of cobases with matching invariants */
		index_set_list matches;

		if ( opts.gramVec ) {

			/* get set of cobases with matching gram vectors */
			cobasis_gram_range range =
					grams.equal_range( fastGramVec(gramMat, cob) );

			/* check invariants for each of these cobases */
			for (cobasis_gram_map::const_iterator it = range.first;
					it != range.second; ++it) {

				vertex_data_ptr oldDat = it->second.second;

				if ( dat->inc.count() == oldDat->inc.count()
						&& dat->gram == oldDat->gram ) {
					matches.push_back(it->second.first);
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

	dfs::vertex_data_list dfs::explorer::matchingInvariants(
			dfs::coordinates_map verts, dfs::vertex_gram_map grams,
			dfs::vertex_data_ptr rep) {

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
		
		initialCobasis = cob->cob;
		ex.cobasisCache.insert(initialCobasis);

		#pragma omp master
		{
		if ( opts.printTrace ) {
			opts.output() << "#I initial basis: " << fmt( cob->cob ) 
					<< " " << *sol << "\n";
		}

		addVertex(dat);
		getRays(ex);
		} /* omp master */
		
		/* DFS the edge graph */
		//TODO FIXME res = dfsFromRoot(cob);
		pushNewEdges(ex, cob->cob);


		nSuccess += true;
		} /* omp parallel */
		
		/* set algorithm end time */
		diff_time = std::clock() - start_time;
		
		return ( nThreads == nSuccess );
	}
	
	////////////////////////////////////////////////////////////////////
	// Query methods for after completion of doDfs()
	////////////////////////////////////////////////////////////////////
	
	dfs::cobasis_map dfs::getBasisOrbits() const {
		return cobasis_map(globalBasisOrbits.begin(), 
				globalBasisOrbits.end());
	}
	
	ind dfs::getDimension() const { return dim - 1; }
	
	index_set dfs::getInitialCobasis() const { return initialCobasis; }
	
	bool dfs::isFinished() const { return !hitMaxBasis; }
	
	dfs::coordinates_map dfs::getRayOrbits() const { 
		return coordinates_map(globalRayOrbits.begin(), 
				globalRayOrbits.end());
	}
	
	std::clock_t dfs::getRunningTime() const 
		{ return diff_time / clocks_per_ms; }
	
	permutation_group const& dfs::getSymmetryGroup() const { return g; }
	
	dfs::coordinates_map dfs::getVertexOrbits() const { 
		return coordinates_map(globalVertexOrbits.begin(), 
				globalVertexOrbits.end());
	}
	
	gram_matrix const& dfs::getGramMat() const { return gramMat; }
	
	
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
		globalBasisOrbits.push_back(std::make_pair(cob, dat));
		oSize = globalBasisOrbits.size();
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
		globalVertexOrbits.push_back(std::make_pair(dat->coords, dat));
		oSize = globalVertexOrbits.size();
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
	
	void dfs::getRays(dfs::explorer& ex) {
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
					if ( ex.rayUpdate == globalRayOrbits.size() ) {
						/* local up to date, add to global */
						globalRayOrbits.push_back(
								std::make_pair(dat->coords, dat));
					} else {
						/* update local to current global state */
						newRayOrbits.insert(
								globalRayOrbits.begin() + ex.rayUpdate,
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
						++ex.rayUpdate;

						if ( opts.printRay && oSize % opts.printRay == 0 ) {
							#pragma omp critical(print)
							{
							std::ostream& out = opts.output();
							out << "# rays: " << oSize << " ("
									<< currentTime() << " ms)";
							if ( opts.printNew ) {
								out << " " << dat->coords;
								if ( opts.debugGram ) out << " " << dat->gram;
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
						ex.rayUpdate += newRayOrbits.size();

						/* Not a new ray, by new local cache */
						if ( ex.knownRay(newRayOrbits, dat) ) break;
					}
				}
			}
		}
	}
	
	void dfs::initGlobals() {
		/* account for flipability of arrangement gram matrix */
		if ( opts.aRepresentation ) gramMat = gramMat.abs();
		
		/* Default initialize remaining data members */
		globalBasisOrbits = cobasis_list();
		hitMaxBasis = false;
		initialCobasis = index_set();
		globalRayOrbits = coordinates_list();
		realDim = 0;
		globalVertexOrbits = coordinates_list();
		globalWorkStack = std::deque<pivot>();
	}
	
	bool dfs::knownOrAddNewCobasis(dfs::explorer& ex,
			index_set cob, dfs::vertex_data_ptr dat) {

		/* check cobasis against local store */
		bool known =
				ex.isKnownCobasis(ex.basisOrbits, ex.cobasisGramMap, cob, dat);
		if ( known ) { return false; }

		while ( true ) {
			cobasis_map newBasisOrbits;
			cobasis_gram_map newCobasisGramMap;

			#pragma omp critical(cobases)
			{
			if ( ex.cobasisUpdate == globalBasisOrbits.size() ) {
				/* local up to date, add to global */
				addCobasis(cob, dat);
			} else {
				/* update local to current global state */
				newBasisOrbits.insert(
						globalBasisOrbits.begin() + ex.cobasisUpdate,
						globalBasisOrbits.end());
			}
			} /* omp critical(cobases) */

			if ( newBasisOrbits.empty() ) {
				/* found a new cobasis */

				vertex_data_ptr newDat = boost::make_shared<vertex_data>(
						dat->coords, dat->inc, dat->cobs, dat->det, dat->gram);

				ex.basisOrbits.insert(std::make_pair(cob, newDat));
				if ( opts.gramVec ) {
					ex.cobasisGramMap.insert(std::make_pair(
							fastGramVec(gramMat, cob),
							std::make_pair(cob, newDat)));
				}
				++ex.cobasisUpdate;

				return true;
			} else {
				/* update locals with new globals */

				for (cobasis_map::iterator it = newBasisOrbits.begin();
						it != newBasisOrbits.end(); ++it) {
					vertex_data& val = *(it->second);
					vertex_data_ptr newDat = boost::make_shared<vertex_data>(
							val.coords, val.inc, val.cobs, val.det, val.gram);

					ex.basisOrbits.insert(std::make_pair(cob, newDat));
					if ( opts.gramVec ) {
						gram_matrix gram = fastGramVec(gramMat, cob);
						ex.cobasisGramMap.insert(
								std::make_pair(gram,
										std::make_pair(cob, newDat)));
						newCobasisGramMap.insert(
								std::make_pair(gram,
										std::make_pair(cob, newDat)));
					}
				}
				ex.cobasisUpdate += newBasisOrbits.size();

				/* test vertex against new local cache */
				known = ex.isKnownCobasis(
						newBasisOrbits, newCobasisGramMap, cob, dat);
				if ( known ) { return false; }
			}
		}
	}

	dfs::vertex_data_known dfs::knownOrAddNewVertex(dfs::explorer& ex,
			dfs::vertex_data_ptr rep) {

		/* check vertex against local store */
		vertex_data_ptr known =
				ex.knownVertex(ex.vertexOrbits, ex.vertexGramMap, rep);
		if ( known ) { return vertex_data_known(known, false); }

		while ( true ) {
			coordinates_map newVertexOrbits;
			vertex_gram_map newVertexGramMap;
			uind oSize;

			#pragma omp critical(vertices)
			{
			if ( ex.vertexUpdate == globalVertexOrbits.size() ) {
				/* local up to date, add to global */
				addVertex(rep);
			} else {
				/* update local to current global state */
				newVertexOrbits.insert(
						globalVertexOrbits.begin() + ex.vertexUpdate,
						globalVertexOrbits.end());
			}
			oSize = globalVertexOrbits.size();
			} /* omp critical(vertices) */

			if ( newVertexOrbits.empty() ) {
				/* found a new vertex */

				vertex_data_ptr dat = boost::make_shared<vertex_data>(
						rep->coords, rep->inc, rep->cobs, rep->det, rep->gram);
				ex.vertexOrbits.insert(std::make_pair(dat->coords, dat));
				if ( opts.gramVec ) {
					ex.vertexGramMap.insert(std::make_pair(dat->gram, dat));
				}
				++ex.vertexUpdate;

				return vertex_data_known(dat, true);
			} else {
				/* update locals with new globals */

				for (coordinates_map::iterator it = newVertexOrbits.begin();
						it != newVertexOrbits.end(); ++it) {
					vertex_data& val = *(it->second);
					vertex_data_ptr dat = boost::make_shared<vertex_data>(
							val.coords, val.inc, val.cobs, val.det, val.gram);

					ex.vertexOrbits.insert(std::make_pair(dat->coords, dat));
					if ( opts.gramVec ) {
						ex.vertexGramMap.insert(std::make_pair(dat->gram, dat));
						newVertexGramMap.insert(std::make_pair(dat->gram, dat));
					}
				}
				ex.vertexUpdate += newVertexOrbits.size();

				/* test vertex against new local cache */
				known = ex.knownVertex(newVertexOrbits, newVertexGramMap, rep);
				if ( known ) { return vertex_data_known(known, false); }
			}
		}
	}

	void dfs::pushNewEdges(explorer& ex, index_set& oldCob) {

		/* for each index in the old cobasis */
		for (index_set_iter it = lrs::begin(oldCob);
				it != lrs::end(oldCob); ++it) {

			/* the leaving index */
			ind leave = *it;
			/* the appropriate entering indices */
			index_set entering(oldCob.size());
			/* the entering index */
			ind enter;

			if ( opts.aRepresentation ) {
				/* use arrangement pivot selection */
				entering = ex.l.arrangementRatio(leave);
			} else if ( opts.lexOnly ) {
				/* calculate entering index lexicographically (BAD) */
				enter = ex.l.lexRatio(leave);
				if (enter >= 0) entering.set(enter); else continue;
			} else {
				/* calculate set of valid entering indices */
				entering = ex.l.allRatio(leave);
			}

			if ( opts.printTrace ) {
				#pragma omp critical(print)
				{
				opts.output() << "#I for leaving index { " << leave
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
				if ( opts.showsAllDicts ) {
					#pragma omp critical(print)
					{
					opts.output() << "\nPivot: " << leave << "=>" << enter;
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

						#pragma omp critical(stacks)
						{
						globalWorkStack.push_back(pivot(oldCob, leave, enter));
						} /* omp critical(stacks) */

						if ( opts.printTrace ) {
							#pragma omp critical(print)
							{
							opts.output() << "#I pushing new vertex: "
									<< fmt( cob->cob ) << " " << *sol << "\n";
							} /* omp critical(print) */
						}
					} else if ( dat->coords == vert.first->coords
							|| ! opts.dualFacetTrick ) {

						/* if this is a new cobasis for a previously seen
						 * vertex, and we are not employing the dual facet
						 * trick to prune the search tree, add the cobasis to
						 * the search stack if it is unique */

						if ( knownOrAddNewCobasis(ex, cob->cob, vert.first) ) {
							#pragma omp critical(stacks)
							{
							globalWorkStack.push_back(
									pivot(oldCob, leave, enter));
							} /* omp critical(stacks) */

							if ( opts.printTrace ) {
								#pragma omp critical(print)
								{
								opts.output() << "#I pushing new cobasis: "
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
						if ( opts.printTrace ) {
							#pragma omp critical(print)
							{
							opts.output() << "#I ignoring cobasis "
									<< fmt( cob->cob ) << " by dual facet "
									"trick\n";
							} /* omp critical(print) */
						}
					}
				} else if ( opts.printTrace ) {
					#pragma omp critical(print)
					{
					opts.output() << "#I seen cobasis " << fmt( cob->cob )
							<< " before\n";
					} /* omp critical(print) */
				}
			}
		}
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
				fastGramVec(gramMat, inc) : gram_matrix();
		
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
