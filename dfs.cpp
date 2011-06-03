#include <algorithm>

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
	
	void dfs::addVertex(dfs::vertex_representation_ptr rep) {
		/* TODO write me */
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
		l.getFirstBasis();
		
		if (opts.showAllDicts) l.printDict();
		
		cobasis_ptr cob(l.getCobasis(0));
		coordinates_ptr sol(l.getVertex());
		cobasis_invariants_ptr rep(cobasisInvariants(cob, sol));
		
		initialCobasis = rep;
		addCobasis(rep);
		addVertex(vertexRep(cob, sol));
		
		/* TODO finish me */
		
	}
	
	void dfs::initGlobals() {
		basisCount = 0;
	}
	
	dfs::vertex_representation_ptr dfs::vertexRep(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		/* TODO lots of stuf in dfs.gap VertexRep() that could be added */
		
		//concatenate the cobasis and extra incidence of the cobasis invariants
		index_list inc(cob->cob.size() + cob->extraInc.size());
		inc.insert(inc.end(), cob->cob.begin(), cob->cob.end());
		inc.insert(inc.end(), cob->extraInc.begin(), cob->extraInc.end());
		//and sort the resulting list
		std::sort(inc.begin(), inc.end());
		
		vertex_representation_ptr rep(
			new vertex_representation(inc, *coords, cob->det)
		);
		
		return rep;
	}
	
} /* namespace basil */
