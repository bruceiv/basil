#include <algorithm>
#include <set>
#include <vector>

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
		
		dfsFirstBasis();
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
		vertexSet.insert(rep->coords / rep->coords[1]);
	}
	
	dfs::cobasis_invariants_ptr dfs::cobasisInvariants(dfs::cobasis_ptr cob, 
			coordinates_ptr coords) { 
		
		cobasis_invariants_ptr cobI(
				new cobasis_invariants(cob->cob, cob->extraInc, *coords, 
									   cob->det)
				);
		
		return cobI;
	}
	
	void dfs::dfsFirstBasis() {
		
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
		rayOrbits = std::vector<vertex_rep_ptr>();
		vertexOrbits = std::vector<vertex_rep_ptr>();
		vertexSet = std::set<coordinates>();
	}
	
	dfs::vertex_rep_ptr dfs::knownRay(dfs::vertex_rep_ptr rep) {
		//TODO write me
		return rep;
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
