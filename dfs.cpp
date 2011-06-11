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
	
	dfs::results dfs::doDfs() {
		initGlobals();
		
		if (opts.showsAllDicts) l.printDict();
		
		index_set cob = dfsFirstBasis()->cob;
		
		l.setCobasis(cob);
		
		bool finished = dfsFromRoot(cob);
		
		return results(cobasisList, dim-1, initialCobasis->cob, finished, 
					   rayOrbits, g, vertexOrbits);
	}
	
	/** Prints a representation of its cobasis (as a set of indices). */
	void print(std::ostream& o, lrs::index_set const& s) {
		bool isFirst = true;
		o << "{";
		for (lrs::index_set_iter it = lrs::begin(s); 
			 it != lrs::end(s); 
			 ++it) {
			if (isFirst) isFirst = false; else o << ", ";
			o << *it;
		}
		o << "}";
	}
	
	void print(std::ostream& o, permutation_list& l) {
		bool isFirst = true;
		o << "{";
		for (permutation_list::const_iterator it = l.begin();
				it != l.end(); ++it) {
			if (isFirst) isFirst = false; else o << ", ";
			o << **it;
		}
		o << "}";
	}
	
	std::ostream& operator<< (std::ostream& o, dfs::results& r) {
		o << "{\n\tdimension: " << r.dimension
				<< "\n\tinitial cobasis: ";
		print(o, r.initialCobasis);
		o << "\n\tsymmetry generators: ";
		print(o, r.symmetryGroup.S);
		o << "\n\tbasis orbits #: " << r.basisOrbits.size()
				<< "\n\tvertex orbits #: " << r.vertexOrbits.size()
				<< "\n\tray orbits #: " << r.rayOrbits.size()
				<< "\n}" << std::endl;
		
		return o;
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
	
	bool dfs::dfsFromRoot(dfs::index_set& root) {
		
		pushNewEdges(root);
		
		while ( ! workStack.empty() && opts.basisLimit > basisCount ) {
			
			cobasis_ptr dict(l.getCobasis(0));
			
			pivot p = workStack.back();
			workStack.pop_back();
			
			while ( dict->cob != p.cob && ! pathStack.empty() ) {
				
				pivot btPivot = pathStack.back();
				pathStack.pop_back();
				l.pivot(btPivot.leave, btPivot.enter);
				dict.reset(l.getCobasis(0));
			}
			
			l.pivot(p.leave, p.enter);
			
			cobasis_ptr cob(l.getCobasis(0));
			getRays();
			pushNewEdges(cob->cob);
			
			pathStack.push_back(
					pivot(index_set(root.size()), p.enter, p.leave) );
		}
		
		/* Did this finish, or terminate at too many bases? */
		return opts.basisLimit > basisCount;
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
		}
		
		return false;
	}

	void dfs::initGlobals() {
		allIndices = index_set(rows+1).set(); /* set all bits of this index set */
		cobasisCache.resize(opts.cacheSize); /* resize the cobasis cache to its 
											  * proper value */
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
		for (vertex_rep_list::iterator it = rayOrbits.begin();
				it != rayOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = (*it)->inc;
			
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
		
		coordinates norm = rep->coords.normalization();
		if ( vertexSet.find(norm) != vertexSet.end() ) {
			/* duplicate vertex */
			return rep;
		}
		
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
		
		cobasis_invariants_list matches;
		
		for (cobasis_invariants_list::iterator it = cobasisList.begin(); 
				it != cobasisList.end(); ++it) {
			
			cobasis_invariants_ptr old = *it;
			
			if ( old->det != rep->det ) continue;
			
			/* if we reach here, all invariant checks have passed. */
			matches.push_back(old);
		}
		
		return matches;
	}

	void dfs::pushNewEdges(dfs::index_set& oldCob) {
		/* TODO add capability for turning off lexOnly option */
		
		for (index_set_iter it = lrs::begin(oldCob); 
				it != lrs::end(oldCob); ++it) {
			
			ind leave = *it;
			index_set entering(oldCob.size());
			
			ind enter = l.lexRatio(leave);
			if (enter >= 0) entering.set(enter);
			
			for (index_set_iter it2 = lrs::begin(entering); 
					it2 != lrs::end(entering); ++it2) {
				
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
