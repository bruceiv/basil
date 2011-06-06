#include <algorithm>
#include <deque>
#include <set>
#include <vector>

#include <permlib/permlib_api.h>

#include "basilCommon.hpp"
#include "dfs.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

namespace basil {
	
	////////////////////////////////////////////////////////////////////////////
	// Public Members
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::doDfs() {
		initGlobals();
		
		if (opts.showAllDicts) l.printDict();
		
		index_list cob = dfsFirstBasis()->cob;
		
		l.setCobasis(cob);
		
		dfsFromRoot(cob);
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Private Members
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::addCobasis(dfs::cobasis_invariants_ptr cob) {
		/* TODO lots of stuff in dfs.gap AddCobasis() that could be added */
		
		cob->index = ++basisCount;
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
		if (opts.showAllDicts) l.printDict();
		
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
	
	void dfs::initGlobals() {
		basisCount = 0;
		cobasisQueue = std::deque<index_list>();
		rayOrbits = std::vector<vertex_rep_ptr>();
		vertexOrbits = std::vector<vertex_rep_ptr>();
		vertexSet = std::set<coordinates>();
	}
	
	dfs::vertex_rep_ptr dfs::knownRay(dfs::vertex_rep_ptr rep) {
		
		/* incidence set to find */
		index_list& find = rep->inc;
		
		/* for every known orbit representative */
		for (std::vector<vertex_rep_ptr>::iterator it = rayOrbits.begin();
				it != rayOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_list& old = (*it)->inc; 
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the ray we are trying to find to the incidence 
			 * set of the known ray. */
			permutation_ptr act = permlib::setImage(
				g, find.begin(), find.end(), old.begin(), old.end());
			
			/* if such a permuation is found, return the known ray */
			if (act) return *it;
		}
		
		/* no known ray that is equivalent up to symmetry */
		return vertex_rep_ptr();
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
				if (opts.showAllDicts) l.printDict();
				cobasis_ptr cob(l.getCobasis(0));
				l.pivot(enter, leave);
				
				
				
				//TODO finish me
			}
		}
	}
	
	dfs::vertex_rep_ptr dfs::rayRep(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		
		/* TODO add gramVec option */
		
		//concatenate the cobasis and extra incidence of the cobasis invariants
		index_list inc(cob->cob.size() + cob->extraInc.size());
		inc.insert(inc.end(), cob->cob.begin(), cob->cob.end());
		inc.insert(inc.end(), cob->extraInc.begin(), cob->extraInc.end());
		//remove the ray index
		inc.erase(std::remove(inc.begin(), inc.end(), cob->ray));
		//and sort the resulting list
		std::sort(inc.begin(), inc.end());
		
		vertex_rep_ptr rep(
			new vertex_rep(inc, *coords, cob->det)
		);
		
		return rep;
	}
	
	dfs::vertex_rep_ptr dfs::vertexRep(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		/* TODO lots of stuf in dfs.gap VertexRep() that could be added */
		
		//concatenate the cobasis and extra incidence of the cobasis invariants
		index_list inc(cob->cob.size() + cob->extraInc.size());
		inc.insert(inc.end(), cob->cob.begin(), cob->cob.end());
		inc.insert(inc.end(), cob->extraInc.begin(), cob->extraInc.end());
		//and sort the resulting list
		std::sort(inc.begin(), inc.end());
		
		vertex_rep_ptr rep(
			new vertex_rep(inc, *coords, cob->det)
		);
		
		return rep;
	}
	
} /* namespace basil */
