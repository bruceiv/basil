#include <algorithm>
#include <ctime>
#include <deque>
#include <iterator>
#include <ostream>
#include <set>
#include <vector>

#include <boost/make_shared.hpp>

#include <permlib/permlib_api.h>
#include <permlib/permutation.h>

#include "dfs.hpp"
#include "fmt.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"

/* simple flag to find debug code in source */
#define BAS_DEBUG 
//#define BAS_DEBUG //

namespace basil {
	
	////////////////////////////////////////////////////////////////////////////
	//
	//  Public Members
	//
	////////////////////////////////////////////////////////////////////////////
	
	bool dfs::doDfs() {
		
		/* set algorithm start time */
		start_time = std::clock();
		
		/* print initial dictionary */
		if (opts.showsAllDicts) l.printDict();
		
		/* get initial cobasis / vertex either from options or from LRS */
		index_set cob = dfsFirstBasis();
		
		/* DFS the edge graph, returning whether it successfully completes */
		bool res = dfsFromRoot(cob);
		
		/* set algorithm end time */
		diff_time = std::clock() - start_time;
		
		return res;
	}

	////////////////////////////////////////////////////////////////////////////
	// Query methods for after completion of doDfs()
	////////////////////////////////////////////////////////////////////////////
	
	dfs::cobasis_map const& dfs::getBasisOrbits() const 
		{ return basisOrbits; }
	
	ind dfs::getDimension() const 
		{ return dim - 1; }
	
	dfs::index_set dfs::getInitialCobasis() const 
		{ return initialCobasis; }
	
	bool dfs::isFinished() const 
		{ return !hitMaxBasis; }
	
	dfs::coordinates_map const& dfs::getRayOrbits() const 
		{ return rayOrbits; }
	
	std::clock_t dfs::getRunningTime() const 
		{ return diff_time / clocks_per_ms; }
	
	permutation_group const& dfs::getSymmetryGroup() const 
		{ return g; }
	
	dfs::coordinates_map const& dfs::getVertexOrbits() const 
		{ return vertexOrbits; }
	
	
	////////////////////////////////////////////////////////////////////////////
	//
	//  Private Members
	//
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::addCobasis(dfs::index_set const& cob, dfs::vertex_data_ptr dat) {
		/* TODO lots of stuff in dfs.gap AddCobasis() that should be looked at 
		 */
		
		basisOrbits.insert(std::make_pair(cob, dat));
		
		/* print cobasis, if option set */
		if ( opts.printBasis && basisOrbits.size() % opts.printBasis == 0 ) {
			std::ostream& out = opts.output();
			out << "# cobases: " << basisOrbits.size() << " (" 
					<< currentTime() << " ms)";
			if ( opts.printNew ) out << " " << fmt( cob );
			BAS_DEBUG out << "\n%" << cob << "%";
			out << std::endl;
		}
	}
	
	void dfs::addVertex(dfs::vertex_data_ptr dat) {
		/* TODO look into handling invariants, gramVec in this framework.
		 * dfs.gap's AddCobasis() and AddVertex() would be helpful */
		
		/* map the normalization of the coordinates to the vertex data */
		coordinates norm = dat->coords.normalization();
		vertexOrbits.insert(std::make_pair(norm, dat));
		
		/* for each defined cobasis, map it to the vertex data */
		for (std::set<index_set>::iterator it = dat->cobs.begin();
				it != dat->cobs.end(); ++it) {
			addCobasis(*it, dat);
		}
		
		/* print vertex, if option set */
		if ( opts.printVertex && vertexOrbits.size() % opts.printVertex == 0 ) {
			std::ostream& out = opts.output();
			out << "# vertices: " << vertexOrbits.size() << " (" 
					<< currentTime() << " ms)";
			if ( opts.printNew ) out << " " << norm;
			BAS_DEBUG out << "\n^" << dat->inc << "^";
			out << std::endl;
		}
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
		vertex_data_ptr dat(vertexData(cob, sol));
		
		initialCobasis = cob->cob;
		addVertex(dat);
		getRays();
		
		/* return the initial cobasic indices */
		return initialCobasis;
	}
	
	bool dfs::dfsFromRoot(basil::dfs::index_set& root) {
		
		/* Add new vertex representations / cobases adjacent to the root 
		 * vertex to the work stack */
		pushNewEdges(root);
		
		/* While there are new (up to symmetry) vertices to explore, and the 
		 * maximum problem size has not been exceeded... */
		while ( ! workStack.empty() && opts.basisLimit > basisOrbits.size() ) {
			
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
		hitMaxBasis = opts.basisLimit <= basisOrbits.size();
		return ! hitMaxBasis;
	}
	
	void dfs::getRays() {
		for (ind j = 1; j <= realDim; j++) {
			coordinates_ptr s( l.getSolution(j) );
			
			if (s) {
				cobasis_ptr c( l.getCobasis(j) );
				vertex_data_ptr dat( rayData(c, s) );
				if ( ! knownRay(dat) ) {
					rayOrbits.insert(std::make_pair(dat->coords, dat));
					
					if ( opts.printRay 
							&& rayOrbits.size() % opts.printRay == 0 ) {
						std::ostream& out = opts.output();
						out << "# rays: " << rayOrbits.size() << " (" 
								<< currentTime() << " ms)";
						if ( opts.printNew ) out << " " << dat->coords;
						BAS_DEBUG out << "\n*" << dat->inc << "*";
						BAS_DEBUG for (std::set<index_set>::iterator it 
						BAS_DEBUG 		= dat->cobs.begin(); 
						BAS_DEBUG		it != dat->cobs.end(); ++it) {
						BAS_DEBUG 	out << "\n%" << *it << "%";
						BAS_DEBUG }
						out << std::endl;
					}
				}
			}
		}
	}
	
	bool dfs::findSymmetry(dfs::index_set find, dfs::index_set_list list) {
		
		/* for each possible size of superset of this cobasis */
		for (ind groundSize = find.count() + 1; groundSize <= rows; 
				groundSize++) {
			
			/* for each cobasis in the list to check for symmetry */
			for (index_set_list::iterator it = list.begin(); it != list.end(); 
					++it) {
				
				/* the cobasis to check for symmetry */
				index_set old = *it;
				
				if ( find == old ) {
					/* duplicate cobasis */
					return true;
				}
				
				if ( opts.assumesNoSymmetry ) continue;
				
				/* Take the set union of the two cobases into ground */
				index_set ground = find | old;
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
				permutation_ptr act = permlib::setImage(*stab, 
						plBegin(find), plEnd(find), plBegin(old), plEnd(old));
				
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
		
		basisOrbits = cobasis_map();
		cobasisQueue = std::deque<index_set>();
		hitMaxBasis = false;
		initialCobasis = index_set();
		pathStack = std::deque<index_pair>();
		rayOrbits = coordinates_map();
		realDim = 0;
		vertexOrbits = coordinates_map();
		workStack = std::deque<pivot>();
		
	}
	
	bool dfs::isNewCobasis(dfs::index_set cob, dfs::vertex_data_ptr dat) {
		
		/* TODO look at adding canonTest, gramMotion from dfs.gap 
		 * IsNewCobasis() */
		
		index_set_list possibleMatches = matchingInvariants(dat);
		
		/* if no known cobasis has invariants matching this one, it's new */
		if ( possibleMatches.size() == 0 ) return true;
		
		/* if a known cobasis (with matching invariants) is symmetric to this 
		 * one, it's not new */
		if ( findSymmetry(cob, possibleMatches) ) return false;
		
		/* if we can't find a matching cobasis, this one must be new */
		return true;
	}

	dfs::vertex_data_ptr dfs::knownRay(dfs::vertex_data_ptr rep) {
		
		/* incidence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (coordinates_map::iterator it = rayOrbits.begin(); 
				it != rayOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = it->second->inc;
			
			/* if we assume no symmetry, simply check for equal cobases */
			if ( opts.assumesNoSymmetry ) {
				if ( find == old ) return it->second; else continue;
			}
			
			/* PermLib chokes on trying to find set image of sets with 
			 * non-equal size. As equal size of incidence sets should be a 
			 * cheap invariant, though, I'll check it here. */
			if ( find.count() != old.count() ) continue;
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the ray we are trying to find to the incidence 
			 * set of the known ray. */
			permutation_ptr act = permlib::setImage(
					g, plBegin(find), plEnd(find), plBegin(old), plEnd(old));
			
			/* if such a permuation is found, return the known ray */
			if (act) return it->second;
		}
		
		/* no known ray that is equivalent up to symmetry */
		return vertex_data_ptr();
	}
	
	dfs::vertex_data_ptr dfs::knownVertex(dfs::vertex_data_ptr rep) {
		
		/* TODO add gramVec / invariants handling */
		
		/* a normalization of this vertex */
		coordinates norm = rep->coords.normalization();
		/* if it's already in the set, it's not new */
		if ( vertexOrbits.find(norm) != vertexOrbits.end() ) {
			/* duplicate vertex */
			return rep;
		}
		
		/* if we assume no symmetry, it must be new */
		if ( opts.assumesNoSymmetry ) return vertex_data_ptr();
		
		/* incedence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
		for (coordinates_map::iterator it = vertexOrbits.begin(); 
				it != vertexOrbits.end(); ++it) {
			
			/* incidence set to check */
			index_set& old = it->second->inc;
			
			/* PermLib chokes on trying to find set image of sets with 
			 * non-equal size. As equal size of incidence sets should be a 
			 * cheap invariant, though, I'll check it here. */
			if ( find.count() != old.count() ) continue;
			
			/* look for a permutation in the global group that maps the 
			 * incidence set of the vertex we are trying to find to the 
			 * incidence set of the known vertex. */
			permutation_ptr act = permlib::setImage(
					g, plBegin(find), plEnd(find), plBegin(old), plEnd(old));
			
			/* if such a permuation is found, return the known vertex */
			if (act) return it->second;
		}
		
		/* no known vertex that is equivalent up to symmetry */
		return vertex_data_ptr();
	}
	
	dfs::index_set_list dfs::matchingInvariants(dfs::vertex_data_ptr dat) {
		
		/* TODO lots of stuff in equivalent Symbal code to add */
		
		/* list of cobases with matching invariants */
		index_set_list matches;
		
		/* for each known cobasis */
		for (cobasis_map::iterator it = basisOrbits.begin();
				it != basisOrbits.end(); ++it) {
			
			vertex_data_ptr old = it->second;
			
			/* skip if they have non-matching determinants */
			if ( old->det != dat->det ) continue;
			
			/* if we reach here, all invariant checks have passed, add the 
			 * cobasis to the list, then */
			matches.push_back(it->first);
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
				 * seen cobases (which could be reached from different pivots).
				 * 
				 * As the time to insert or lookup should be similar, using 
				 * insert instead of lookup saves the extra call to insert here 
				 */
				if ( ! cobasisCache.insert(cob->cob) ) {
					
					/* if this cobasis is not in the cache, add it */
					/* cobasisCache.insert(cob->cob); */
					
					/* calculate invariants of new cobasis */
					vertex_data_ptr newVertex(vertexData(cob, sol));
					vertex_data_ptr oldVertex(knownVertex(newVertex));
					
					if ( ! oldVertex ) {
						
						/* this vertex has yet to be seen, add it */
						addVertex(newVertex);
						/* add this vertex to the search stack */
						workStack.push_back(pivot(oldCob, leave, enter));
						
					} else if (oldVertex->coords == newVertex->coords 
								|| ! opts.dualFacetTrick ) {
						
						/* if this is a new cobasis for a previously seen 
						 * vertex, and we are not employing the dual facet 
						 * trick to prune the search tree, add the cobasis to 
						 * the search stack if it is unique */
						if ( isNewCobasis(cob->cob, newVertex) ) {
							addCobasis(cob->cob, oldVertex);
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
	
	dfs::vertex_data_ptr dfs::rayData(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		/* TODO look at including gramVec from dfs.gap RayRep() */
		
		/* union of the cobasis and extra incidence of the cobasis data */
		index_set inc = cob->cob | cob->extraInc;
		/* less the ray index */
		inc.set(cob->ray, false);
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				*coords, inc, cob->cob, cob->det);
		
		return dat;
	}
	
	dfs::vertex_data_ptr dfs::vertexData(dfs::cobasis_ptr cob, 
			dfs::coordinates_ptr coords) {
		/* TODO look at including gramVec, stabilizerOrbits from dfs.gap 
		 * VertexRep() */
		
		/* union of the cobasis and extra incidence of the cobasis data */
		index_set inc = cob->cob | cob->extraInc;
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				*coords, inc, cob->cob, cob->det);
		
		return dat;
	}
	
} /* namespace basil */
