#include <algorithm>
#include <deque>
#include <iterator>
#include <set>
#include <vector>

#include <permlib/permlib_api.h>

#include "basilCommon.hpp"
#include "dfs.hpp"
#include "lruCache.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

namespace basil {
	
	////////////////////////////////////////////////////////////////////////////
	// Public Members
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::doDfs() {
		initGlobals();
		
		if (opts.showsAllDicts) l.printDict();
		
		index_set cob = dfsFirstBasis()->cob;
		
		l.setCobasis(cob);
		
		dfsFromRoot(cob);
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Private Members
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::addCobasis(dfs::cobasis_invariants_ptr cob) {
		/* TODO lots of stuff in dfs.gap AddCobasis() that could be added */
		
		cob->index = ++basisCount;
		cobasisList.push_back(cob);
	}
	
	void dfs::addVertex(dfs::vertex_rep_ptr rep) {
		vertexOrbits.push_back(rep);
		
		//TODO consult with Dr. Bremner about the normalization algo
		coordinates norm = rep->coords.normalization();
		vertexSet.insert(norm);
	}
	
	dfs::cobasis_invariants_ptr dfs::cobasisInvariants(dfs::cobasis_ptr cob, 
			coordinates_ptr coords) { 
		
		cobasis_invariants_ptr cobI(
				new cobasis_invariants(cob->cob, cob->extraInc, *coords, 
									   cob->det)
				);
		
		return cobI;
	}
	
	dfs::cobasis_ptr dfs::dfsFirstBasis() {
		
		if (! l.getFirstBasis() ) 
			throw dfs_error("LRS failed to find first basis.");
		if (opts.showsAllDicts) l.printDict();
		
		realDim = l.getRealDim();
		cobasis_ptr cob(l.getCobasis(0));
		coordinates_ptr sol(l.getVertex());
		cobasis_invariants_ptr rep(cobasisInvariants(cob, sol));
		
		initialCobasis = rep;
		addCobasis(rep);
		addVertex(vertexRep(cob, sol));
		getRays();
		
		return cob;
	}
	
	void dfs::dfsFromRoot(dfs::index_list& root) {
		
		pushNewEdges(root);
		
		//TODO finish me
	}
	
	void dfs::getRays() {
		for (ind j = 1; j <= realDim; j++) {
			coordinates_ptr s( l.getSolution(j) );
			
			if (s) {
				cobasis_ptr c( l.getCobasis(j) );
				vertex_rep_ptr rep( rayRep(c, s) );
				if (! knownRay(rep) ) rayOrbits.push_back(rep);
			}
		}
	}
	
	bool dfs::findSymmetry(dfs::cobasis_invariants_ptr rep, 
						   dfs::cobasis_invariants_list list) {
		
		for (ind groundSize = rep->cob.count()+1; groundSize <= rows; 
				groundSize++) {
			
			for (cobasis_invariants_list::iterator it = list.begin(); 
					it != list.end(); ++it) {
				
				cobasis_invariants_ptr old = (*it);
				
				if ( rep->cob == old->cob ) {
					//duplicate cobasis
					return true;
				}
				
				/* Take the set union into ground */
				index_set ground = rep->cob | old->cob;
				/* Take the complement of ground into leftOut */
				index_set leftOut = allIndices - ground;
								
				while ( ground.count() < groundSize ) {
					//take a random left out element and add it to ground
					ind randInd = pseudoRandomInd(leftOut);
					ground.set(randInd, true);
					leftOut.set(randInd, false);
				}
				
				
				
				//TODO finish me
			}
		}
		
		return false;
	}

	void dfs::initGlobals() {
		allIndices = index_set(rows).set(); //an index set with all bits set
		basisCount = 0;
		cobasisCache = lru_cache<index_set>(opts.cacheSize);
		cobasisQueue = std::deque<index_set>();
		cobasisList = cobasis_invariants_list();
		rayOrbits = std::vector<vertex_rep_ptr>();
		vertexOrbits = std::vector<vertex_rep_ptr>();
		vertexSet = std::set<coordinates>();
	}
	
	bool dfs::isNewCobasis(dfs::cobasis_invariants_ptr rep) {
		
		/* TODO lots of options in original symbal that could be added */
		
		cobasis_invariants_list possibleMatches = matchingInvariants(rep);
		
		/* new by invariants */
		if ( possibleMatches.size() == 0 ) return true;
		
		if ( findSymmetry(rep, possibleMatches) ) return false;
		
		return true;
	}

	dfs::vertex_rep_ptr dfs::knownRay(dfs::vertex_rep_ptr rep) {
		
		/* incidence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (std::vector<vertex_rep_ptr>::iterator it = rayOrbits.begin();
				it != rayOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_list& old = (*it)->inc; 
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the ray we are trying to find to the incidence 
			 * set of the known ray. */
			permutation_ptr act = permlib::setImage(
				g, begin(find), end(find), begin(old), end(old));
			
			/* if such a permuation is found, return the known ray */
			if (act) return *it;
		}
		
		/* no known ray that is equivalent up to symmetry */
		return vertex_rep_ptr();
	}
	
	dfs::vertex_rep_ptr dfs::knownVertex(dfs::vertex_rep_ptr rep) {
		
		/* TODO add gramVec / invariants handling */
		
		coordinates norm = rep->coords.normalization();
		if ( vertexSet.find(norm) != vertexSet.end() ) {
			/* duplicate vertex */
			return rep;
		}
		
		/* incedence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (std::vector<vertex_rep_ptr>::iterator it = vertexOrbits.begin(); 
				it != vertexOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = (*it)->inc;
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the vertex we are trying to find to the 
			 * incidence set of the known vertex. */
			permutation_ptr act = permlib::setImage(
				g, begin(find), end(find), begin(old), end(old));
			
			/* if such a permuation is found, return the known ray */
			if (act) return *it;
		}
		
		/* no known vertex that is equivalent up to symmetry */
		return vertex_rep_ptr();
	}
	
	dfs::cobasis_invariants_list dfs::matchingInvariants(
			dfs::cobasis_invariants_ptr rep) {
		
		/* TODO lots of stuff in equivalent Symbal code to add */
		
		cobasis_invariants_list matches;
		
		for (cobasis_invariants_list::iterator it = cobasisList.begin(); 
				it != cobasisList.end(); ++it) {
			
			cobasis_invariants_ptr old = *it;
			
			if ( ! old->det == rep->det ) continue;
			
			/* if we reach here, all invariant checks have passed. */
			matches.push_back(old);
		}
		
		return matches;
	}

	void dfs::pushNewEdges(dfs::index_list& oldCob) {
		/* TODO add capability for turning off lexOnly option */
		
		for (index_list::iterator it = oldCob.begin(); 
				it != oldCob.end(); ++it) {
			
			ind leave = *it;
			index_list entering;
		
			ind enter = l.lexRatio(leave);
			if (enter >= 0) entering.push_back(enter);
			
			for (index_list::iterator it2 = entering.begin(); 
					it2 != entering.end(); ++it) {
				
				enter = *it2;
				
				/* Do the given pivot, then get the cobasis for the new edge */
				l.pivot(leave, enter);
				if (opts.showsAllDicts) l.printDict();
				cobasis_ptr cob(l.getCobasis(0));
				coordinates_ptr sol(l.getVertex());
				l.pivot(enter, leave);
				
				/*TODO verify with Dr. Bremner that just using the cobasis set
				 * of this is acceptable, rather than the whole record */
				if ( ! cobasisCache.lookup(cob->cob) ) {
					
					/* if this cobasis is not in the cache, add it and 
					 * recalculate */
					cobasisCache.insert(cob->cob);
					
					cobasis_invariants_ptr newRep(cobasisInvariants(cob, sol));
					vertex_rep_ptr newVertex(vertexRep(cob, sol));
					vertex_rep_ptr oldVertex(knownVertex(newVertex));
					
					if ( ! oldVertex ) {
						
						/* this vertex has yet to be seen */
						addVertex(newVertex);
						addCobasis(newRep);
						workStack.push_back(pivot(oldCob, leave, enter));
						
					} else if ( oldVertex->coords == newRep->coords  
								|| ! opts.dualFacetTrick ) {
						
						if ( isNewCobasis(newRep) ) {
							addCobasis(newRep);
							workStack.push_back(pivot(oldCob, leave, enter));
						}
						
					} /* else {
						// we assume that if the new cobasis is defining a 
						// different vertex, but that vertex is symmetric, then 
						// its neighbours will be symmetric to those of the 
						// known vertex. Prune via the dual facet trick.
					} */
				}
			}
		}
	}
	
	dfs::vertex_rep_ptr dfs::rayRep(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		
		/* TODO add gramVec option */
		
		//concatenate the cobasis and extra incidence of the cobasis invariants
		index_set inc = cob->cob | cob->extraInc;
		//remove the ray index
		inc.set(cob->ray, false);
		
		vertex_rep_ptr rep(
			new vertex_rep(inc, *coords, cob->det)
		);
		
		return rep;
	}
	
	dfs::vertex_rep_ptr dfs::vertexRep(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		/* TODO lots of stuf in dfs.gap VertexRep() that could be added */
		
		//concatenate the cobasis and extra incidence of the cobasis invariants
		index_set inc = cob->cob | cob->extraInc;
		
		vertex_rep_ptr rep(
			new vertex_rep(inc, *coords, cob->det)
		);
		
		return rep;
	}
	
} /* namespace basil */
