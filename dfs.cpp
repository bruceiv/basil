#include <algorithm>
#include <deque>
#include <iterator>
#include <ostream>
#include <set>
#include <vector>

#include <permlib/permlib_api.h>
#include <permlib/permutation.h>

#include "basilCommon.hpp"
#include "dfs.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"

namespace basil {
	
	////////////////////////////////////////////////////////////////////////////
	// Public Members
	////////////////////////////////////////////////////////////////////////////
	
	bool dfs::doDfs() {
		/* set up algorithm globals */
		initGlobals();
		
		/* print initial dictionary */
		if (opts.showsAllDicts) l.printDict();
		
		/* get initial cobasis / vertex either from options or from LRS */
		index_set cob = dfsFirstBasis();
		
		/* DFS the edge graph, returning whether it successfully completes */
		return dfsFromRoot(cob);
		
// 		return results(cobasisList, dim-1, initialCobasis->cob, !hitMaxBasis, 
// 					   rayOrbits, g, vertexOrbits);
	}

	////////////////////////////////////////////////////////////////////////
	// Query methods for after completion of doDfs()
	////////////////////////////////////////////////////////////////////////
	
	dfs::cobasis_invariants_list const& dfs::getBasisOrbits() const 
		{ return cobasisList; }
	
	ind dfs::getDimension() const 
		{ return dim - 1; }
	
	dfs::index_set dfs::getInitialCobasis() const 
		{ return initialCobasis->cob; }
	
	bool dfs::isFinished() const 
		{ return !hitMaxBasis; }
	
	dfs::vertex_rep_list const& dfs::getRayOrbits() const 
		{ return rayOrbits; }
	
	permutation_group const& dfs::getSymmetryGroup() const 
		{ return g; }
	
	dfs::vertex_rep_list const& dfs::getVertexOrbits() const 
		{ return vertexOrbits; }
	
	
	////////////////////////////////////////////////////////////////////////////
	// Private Members
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::addCobasis(dfs::cobasis_invariants_ptr cob) {
		/* TODO lots of stuff in dfs.gap AddCobasis() that could be added */
		/* TODO see if cob->index is actually ever used ... */
		
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
	
	dfs::index_set dfs::dfsFirstBasis() {
		
		/* get initial solution from LRS */
		if (! l.getFirstBasis() ) 
			throw dfs_error("LRS failed to find first basis.");
		/* print initial solution */
		if (opts.showsAllDicts) l.printDict();
		
		/* go to the initial cobasis, if supplied */
		if ( opts.firstCobasis ) l.setCobasis( *opts.firstCobasis );
		
		/* get true problem dimension */
		realDim = l.getRealDim();
		/* add this vertex / cobasis as a representative of its orbit 
		 * (obviously there aren't any others yet) */
		cobasis_ptr cob(l.getCobasis(0));
		coordinates_ptr sol(l.getVertex());
		cobasis_invariants_ptr rep(cobasisInvariants(cob, sol));
		
		initialCobasis = rep; /* save the initial cobasis */
		addCobasis(rep);
		addVertex(vertexRep(cob, sol));
		getRays();
		
		/* return the initial cobasic indices */
		return cob->cob;
	}
	
	bool dfs::dfsFromRoot(basil::dfs::index_set& root) {
		
		/* Add new vertex representations / cobases adjacent to the root 
		 * vertex to the work stack */
		pushNewEdges(root);
		
		/* While there are new (up to symmetry) vertices to explore, and the 
		 * maximum problem size has not been exceeded... */
		while ( ! workStack.empty() && opts.basisLimit > basisCount ) {
			
			/* get the current cobasis */
			cobasis_ptr dict(l.getCobasis(0));
			
			/* pop the pivot to the edge to explore off the work stack */
			pivot p = workStack.back(); workStack.pop_back();
			
			/* backtrack LRS to a state consistent with the pivot to make */
			while ( dict->cob != p.cob && ! pathStack.empty() ) {
				
				/* get the backtracking pivot off the path stack */
				index_pair btPivot = pathStack.back(); pathStack.pop_back();
				/* reverse the pivot */
				l.pivot(btPivot.second, btPivot.first);
				/* reset the current cobasis */
				dict.reset(l.getCobasis(0));
			}
			
			/* pivot to the vertex to explore */
			l.pivot(p.leave, p.enter);
			/* print the dictionary pivoted to */
			if ( opts.showsAllDicts ) l.printDict();
			
			/* get the new cobasis */
			cobasis_ptr cob(l.getCobasis(0));
			/* get the new rays */
			getRays();
			
			/* Add new vertex representations / cobases adjacent to the new 
			 * vertex to the work stack */
			pushNewEdges(cob->cob);
			
			/* push the pivot just made onto the backtracking stack */
			pathStack.push_back( index_pair(p.leave, p.enter) );
		}
		
		/* Did this finish, or terminate at too many bases? */
		hitMaxBasis = opts.basisLimit <= basisCount;
		return ! hitMaxBasis;
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
		
		/* for each possible size of superset of this cobasis */
		for (ind groundSize = rep->cob.count()+1; groundSize <= rows; 
				groundSize++) {
			
			/* for each cobasis in the list to check for symmetry */
			for (cobasis_invariants_list::iterator it = list.begin(); 
					it != list.end(); ++it) {
				
				/* the cobasis to check for symmetry */
				cobasis_invariants_ptr old = (*it);
				
				if ( rep->cob == old->cob ) {
					/* duplicate cobasis */
					return true;
				}
				
				if ( opts.assumesNoSymmetry ) continue;
				
				/* Take the set union of the two cobases into ground */
				index_set ground = rep->cob | old->cob;
				/* Take the complement of ground into leftOut */
				index_set leftOut = allIndices - ground;
				
				while ( ground.count() < uind(groundSize) ) {
					/* take a random left out element and add it to ground */
					ind randInd = lrs::pseudoRandomInd(leftOut);
					ground.set(randInd, true);
					leftOut.set(randInd, false);
				}
				
				/* get the set stabilizer of ground */
				permutation_group_ptr stab = permlib::setStabilizer(
						g, plBegin(ground), plEnd(ground));
				
				/* look for a permutation in the stabilizer group that maps the 
				 * incidence set of the cobasis we are trying to find to the 
				 * incidence set of the known cobasis. */
				permutation_ptr act = permlib::setImage(
						*stab, plBegin(rep->cob), plEnd(rep->cob), 
						plBegin(old->cob), plEnd(old->cob));
				
				/* This cobasis is symmetric to one we already know of */
				if (act) return true;
			}
			
			/* only check the cobases once if assuming no symmetry */
			if ( opts.assumesNoSymmetry ) break;
		}
		
		/* no symmetry between this cobasis and any other in the list */
		return false;
	}

	void dfs::initGlobals() {
		/* represents the set [1..rows] */
		allIndices = index_set(rows+1).set().set(0, false);
		/* resize the cobasis cache to its proper size */
		cobasisCache.resize(opts.cacheSize);
	}
	
	bool dfs::isNewCobasis(dfs::cobasis_invariants_ptr rep) {
		
		/* TODO lots of options in original symbal that could be added */
		
		cobasis_invariants_list possibleMatches = matchingInvariants(rep);
		
		/* if no known cobasis has invariants matching this one, it's new */
		if ( possibleMatches.size() == 0 ) return true;
		
		/* if a known cobasis (with matching invariants) is symmetric to this 
		 * one, it's not new */
		if ( findSymmetry(rep, possibleMatches) ) return false;
		
		/* if we can't find a matching cobasis, this one must be new */
		return true;
	}

	dfs::vertex_rep_ptr dfs::knownRay(dfs::vertex_rep_ptr rep) {
		
		/* incidence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (vertex_rep_list::iterator it = rayOrbits.begin();
				it != rayOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = (*it)->inc;
			
			/* if we assume no symmetry, simply check for equal cobases */
			if ( opts.assumesNoSymmetry ) {
				if ( find == old ) return *it; else continue;
			}
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the ray we are trying to find to the incidence 
			 * set of the known ray. */
			permutation_ptr act = permlib::setImage(
					g, plBegin(find), plEnd(find), plBegin(old), plEnd(old));
			
			/* if such a permuation is found, return the known ray */
			if (act) return *it;
		}
		
		/* no known ray that is equivalent up to symmetry */
		return vertex_rep_ptr();
	}
	
	dfs::vertex_rep_ptr dfs::knownVertex(dfs::vertex_rep_ptr rep) {
		
		/* TODO add gramVec / invariants handling */
		
		/* a normalization of this vertex */
		coordinates norm = rep->coords.normalization();
		/* if it's already in the set, it's not new */
		if ( vertexSet.find(norm) != vertexSet.end() ) {
			/* duplicate vertex */
			return rep;
		}
		
		/* if we assume no symmetry, it must be new */
		if ( opts.assumesNoSymmetry ) return vertex_rep_ptr();
		
		/* incedence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (vertex_rep_list::iterator it = vertexOrbits.begin(); 
				it != vertexOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = (*it)->inc;
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the vertex we are trying to find to the 
			 * incidence set of the known vertex. */
			permutation_ptr act = permlib::setImage(
					g, plBegin(find), plEnd(find), plBegin(old), plEnd(old));
			
			/* if such a permuation is found, return the known ray */
			if (act) return *it;
		}
		
		/* no known vertex that is equivalent up to symmetry */
		return vertex_rep_ptr();
	}
	
	dfs::cobasis_invariants_list dfs::matchingInvariants(
			dfs::cobasis_invariants_ptr rep) {
		
		/* TODO lots of stuff in equivalent Symbal code to add */
		
		/* list of cobases with matching invariants */
		cobasis_invariants_list matches;
		
		/* for each known cobasis */
		for (cobasis_invariants_list::iterator it = cobasisList.begin(); 
				it != cobasisList.end(); ++it) {
			
			cobasis_invariants_ptr old = *it;
			
			/* skip if they have non-matching determinants */
			if ( old->det != rep->det ) continue;
			
			/* if we reach here, all invariant checks have passed, add the 
			 * cobasis to the list, then */
			matches.push_back(old);
		}
		
		return matches;
	}

	void dfs::pushNewEdges(dfs::index_set& oldCob) {
		
		/* for each index in the old cobasis */
		for (index_set_iter it = lrs::begin(oldCob); 
				it != lrs::end(oldCob); ++it) {
			
			/* the leaving index */
			ind leave = *it;
			/* the appropriate entering indices */
			index_set entering(oldCob.size());
			/* the entering index */
			ind enter;
			
			if (opts.lexOnly) {
				/* calculate entering index lexicographically (BAD) */
				enter = l.lexRatio(leave);
				if (enter >= 0) entering.set(enter); else continue;
			} else {
				/* calculate set of valid entering indices */
				entering = l.allRatio(leave);
			}
			
			/* for each valid entering index */
			for (index_set_iter it2 = lrs::begin(entering); 
					it2 != lrs::end(entering); ++it2) {
				
				enter = *it2;
				
				/* Do the given pivot, then get the cobasis for the new edge */
				l.pivot(leave, enter);
				cobasis_ptr cob(l.getCobasis(0));
				coordinates_ptr sol(l.getVertex());
				/* pivot back */
				l.pivot(enter, leave);
				
				/* avoid expensive invariant calculations by caching recently 
				 * seen cobases (which could be reached from different pivots) 
				 */
				if ( ! cobasisCache.lookup(cob->cob) ) {
					
					/* if this cobasis is not in the cache, add it */
					cobasisCache.insert(cob->cob);
					
					/* calculate invariants of new cobasis */
					cobasis_invariants_ptr newRep(cobasisInvariants(cob, sol));
					vertex_rep_ptr newVertex(vertexRep(cob, sol));
					vertex_rep_ptr oldVertex(knownVertex(newVertex));
					
					if ( ! oldVertex ) {
						
						/* this vertex has yet to be seen, add it */
						addVertex(newVertex);
						addCobasis(newRep);
						/* add this vertex to the search stack */
						workStack.push_back(pivot(oldCob, leave, enter));
						
					} else if ( oldVertex->coords == newRep->coords  
								|| ! opts.dualFacetTrick ) {
						
						/* if this is a new cobasis for a previously seen 
						 * vertex, and we are not employing the dual facet 
						 * trick to prune the search tree, add the cobasis to 
						 * the search stack if it is unique */
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
