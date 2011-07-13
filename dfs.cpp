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
	
	dfs::cobasis_map const& dfs::getBasisOrbits() const { return basisOrbits; }
	
	ind dfs::getDimension() const { return dim - 1; }
	
	dfs::index_set dfs::getInitialCobasis() const { return initialCobasis; }
	
	bool dfs::isFinished() const { return !hitMaxBasis; }
	
	dfs::coordinates_map const& dfs::getRayOrbits() const { return rayOrbits; }
	
	std::clock_t dfs::getRunningTime() const 
		{ return diff_time / clocks_per_ms; }
	
	permutation_group const& dfs::getSymmetryGroup() const { return g; }
	
	dfs::coordinates_map const& dfs::getVertexOrbits() const 
		{ return vertexOrbits; }
	
	basil::matrix const& dfs::getInnerProdMat() const { return innerProdMat; }
	
	dfs::cobasis_gram_map const& dfs::getCobasisGramMap() const 
		{ return cobasisGramMap; }
	
	dfs::vertex_gram_map const& dfs::getVertexGramMap() const 
		{ return vertexGramMap; }

	
	////////////////////////////////////////////////////////////////////////////
	//
	//  Private Members
	//
	////////////////////////////////////////////////////////////////////////////
	
	void dfs::initGlobals() {
		/* represents the set [1..rows] */
		allIndices = index_set(rows+1).set().set(0, false);
		/* resize the cobasis cache to its proper size */
		cobasisCache.resize(opts.cacheSize);
		
		/* Default initialize remaining data members */
		basisOrbits = cobasis_map();
		cobasisGramMap = cobasis_gram_map();
		cobasisQueue = std::deque<index_set>();
		hitMaxBasis = false;
		initialCobasis = index_set();
		pathStack = std::deque<index_pair>();
		rayOrbits = coordinates_map();
		realDim = 0;
		vertexOrbits = coordinates_map();
		vertexGramMap = vertex_gram_map();
		workStack = std::deque<pivot>();
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
		vector_mpz_ptr sol(l.getVertex());
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
				vector_mpz_ptr sol(l.getVertex());
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
	
	void dfs::getRays() {
		for (ind j = 1; j <= realDim; j++) {
			vector_mpz_ptr s( l.getSolution(j) );
			
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
						out << std::endl;
					}
				}
			}
		}
	}
	
	dfs::vertex_data_ptr dfs::knownVertex(dfs::vertex_data_ptr rep) {
		
		/* TODO add gramVec / invariants handling */
		
		/* if it's already in the set, it's not new */
		if ( vertexOrbits.find(rep->coords) != vertexOrbits.end() ) {
			/* duplicate vertex */
			return rep;
		}
		
		/* if we assume no symmetry, it must be new */
		if ( opts.assumesNoSymmetry ) return vertex_data_ptr();
		
		/* List of vertices with matching invariants */
		vertex_data_list possibleMatches = matchingInvariants(rep);
		
		/* New by invariants */
		if ( possibleMatches.size() == 0 ) return vertex_data_ptr();
		
		/* incedence set to find */
		index_set& find = rep->inc;
		
		/* for every known orbit representative */
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
	
	bool dfs::isNewCobasis(dfs::index_set cob, dfs::vertex_data_ptr dat) {
		
		/* TODO look at adding canonTest, gramMotion from dfs.gap 
		 * IsNewCobasis() */
		
		index_set_list possibleMatches = matchingCobasisInvariants(cob, dat);
		
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
			
			/* PermLib throws an assertion failure on setImage of inequally 
			 * sized sets, so I check this very cheap invariant here. */
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
	
	dfs::index_set_list dfs::matchingCobasisInvariants(dfs::index_set cob, 
			dfs::vertex_data_ptr dat) {
		
		/* TODO add cobasis stabilizer invariant from Symbal (?) */
		
		/* list of cobases with matching invariants */
		index_set_list matches;
		
		if (opts.gramVec) {
			
			/* get set of cobases with matching gram vectors */
			cobasis_gram_range range = 
				cobasisGramMap.equal_range( fastGramVec(cob) );
			
			/* check invariants for each of these cobases */
			for (cobasis_gram_map::const_iterator it = range.first; 
					it != range.second; ++it) {
				
				if ( cobasisInvariantsMatch(*it->second.second, *dat) ) {
					matches.push_back(it->second.first);
				}
				
			}
			
		} else {
			
			/* check invariants for each known cobasis */
			for (cobasis_map::const_iterator it = basisOrbits.begin();
					it != basisOrbits.end(); ++it) {
				
				if ( invariantsMatch(*it->second, *dat) ) {
					matches.push_back(it->first);
				}
				
			}
			
		}
		
		return matches;
	}
	
	dfs::vertex_data_list dfs::matchingInvariants(dfs::vertex_data_ptr rep) {
		
		/* list of vertices with matching invariants */
		vertex_data_list matches;
		
		if (opts.gramVec) {
			
			/* get set of cobases with matching gram vectors */
			vertex_gram_range range = 
				vertexGramMap.equal_range( rep->gram );
			
			/* check invariants for each of these cobases */
			for (vertex_gram_map::const_iterator it = range.first; 
					it != range.second; ++it) {
				
				if ( invariantsMatch(*it->second, *rep) ) {
					matches.push_back(it->second);
				}
				
			}
			
		} else {
			
			/* check invariants for each known cobasis */
			for (coordinates_map::const_iterator it = vertexOrbits.begin();
					it != vertexOrbits.end(); ++it) {
				
				if ( invariantsMatch(*it->second, *rep) ) {
					matches.push_back(it->second);
				}
				
			}
			
		}
		
		return matches;
	}
	
	bool dfs::findSymmetry(dfs::index_set find, dfs::index_set_list list) {
		
		if ( opts.assumesNoSymmetry || ! opts.stabSearch ) {
			/* no stabilizer search, so only one pass through the list */
			
			/* for each cobasis in the list to check for symmetry */
			for (index_set_list::iterator it = list.begin(); 
					it != list.end(); ++it) {
				
				/* the cobasis to check for symmetry */
				index_set old = *it;
				
				if ( find == old ) {
					/* duplicate cobasis */
					return true;
				}
				
				/* if no symmetry, check next cobasis */
				if ( opts.assumesNoSymmetry ) continue;
				
				/* look for a permutation in the permutation group that maps 
				 * the incidence set of the cobasis we are trying to find to 
				 * the incidence set of the known cobasis. */
				permutation_ptr act = permlib::setImage(g, 
						plBegin(find), plEnd(find), plBegin(old), plEnd(old));
				
				/* This cobasis is symmetric to one we already know of */
				if (act) return true;
			}
			
		} else {
			/* using stabilizer search, so need to check multiple grounds */
			
			/* for each possible size of superset of this cobasis */
			for (ind groundSize = find.count() + 1; groundSize <= rows; 
					groundSize++) {
				
				/* for each cobasis in the list to check for symmetry */
				for (index_set_list::iterator it = list.begin(); 
						it != list.end(); ++it) {
					
					/* the cobasis to check for symmetry */
					index_set old = *it;
					
					if ( find == old ) {
						/* duplicate cobasis */
						return true;
					}
					
					/* Take the set union of the two cobases into ground */
					index_set ground = find | old;
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
							plBegin(find), plEnd(find), 
							plBegin(old), plEnd(old));
					
					/* This cobasis is symmetric to one we already know of */
					if (act) return true;
				}
			}
		}
		
		/* no symmetry between this cobasis and any other in the list */
		return false;
	}
	
	matrix dfs::fastGramVec(dfs::index_set inc) {
		/* restrict the inner product matrix to the incidence set */
		matrix gram(innerProdMat.restriction(inc));
		
		/* sort each row */
		for (matrix::iterator row = gram.begin(); row != gram.end(); ++row) {
			std::sort( (*row).begin(), (*row).end() );
		}
		
		/* sort the rows */
		std::sort(gram.begin(), gram.end());
		
		return gram;
	}
	
	bool dfs::cobasisInvariantsMatch(dfs::vertex_data const& a, 
			dfs::vertex_data const& b) {
		return true
// 				&& a.det == b.det						/* determinant */
				&& a.inc.count() == b.inc.count()		/* # incident facets */
				&& a.gram == b.gram						/* gram matrix */
				;
	}
	
	bool dfs::invariantsMatch(dfs::vertex_data const& a, 
			dfs::vertex_data const& b) {
		return true
// 				&& a.det == b.det						/* determinant */
				&& a.inc.count() == b.inc.count()		/* # incident facets */
				/* NOTE gram checked by gram matrix lookup for vertices */
				// && a.gram == b.gram					/* gram matrix */
				;
	}
	
	dfs::vertex_data_ptr dfs::rayData(dfs::cobasis_ptr cob, 
			dfs::vector_mpz_ptr coords) {
		/* TODO look at including gramVec from dfs.gap RayRep() */
		
		/* union of the cobasis and extra incidence of the cobasis data */
		index_set inc = cob->cob | cob->extraInc;
		/* less the ray index */
		inc.set(cob->ray, false);
		/* ignore gram vector for the moment */
		matrix gram = matrix(0,0);
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				coordinates(*coords), inc, cob->cob, abs(cob->det), gram);
		
		return dat;
	}
	
	void dfs::addCobasis(dfs::index_set const& cob, dfs::vertex_data_ptr dat) {
		/* TODO lots of stuff in dfs.gap AddCobasis() that should be looked at 
		 */
		
		basisOrbits.insert(std::make_pair(cob, dat));
		
		/* add gram vector, if option set */
		if (opts.gramVec) {
			cobasisGramMap.insert(
				std::make_pair(fastGramVec(cob), std::make_pair(cob, dat))
			);
		}
		
		/* print cobasis, if option set */
		if ( opts.printBasis && basisOrbits.size() % opts.printBasis == 0 ) {
			std::ostream& out = opts.output();
			out << "# cobases: " << basisOrbits.size() << " (" 
					<< currentTime() << " ms)";
			if ( opts.printNew ) out << " " << fmt( cob );
			out << std::endl;
		}
	}
	
	void dfs::addVertex(dfs::vertex_data_ptr dat) {
		/* TODO look into handling invariants, gramVec in this framework.
		 * dfs.gap's AddCobasis() and AddVertex() would be helpful */
		
		/* map the rationalization of the coordinates to the vertex data */
		vertexOrbits.insert(std::make_pair(dat->coords, dat));
		
		/* add gram vector, if option set */
		if (opts.gramVec) {
			vertexGramMap.insert(std::make_pair(fastGramVec(dat->inc), dat));
		}
		
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
			if ( opts.printNew ) out << " " << dat->coords;
			out << std::endl;
		}
	}
	
	dfs::vertex_data_ptr dfs::vertexData(dfs::cobasis_ptr cob, 
			dfs::vector_mpz_ptr coords) {
		/* TODO look at including stabilizerOrbits from dfs.gap VertexRep() */
		
		/* union of the cobasis and extra incidence of the cobasis data */
		index_set inc = cob->cob | cob->extraInc;
		/* gram matrix, or empty if option off */
		matrix gram = ( opts.gramVec ) ? fastGramVec(inc) : matrix(0,0);
		
		vertex_data_ptr dat = boost::make_shared<vertex_data>(
				coords->rationalization(), inc, cob->cob, abs(cob->det), gram);
		
		return dat;
	}
	
} /* namespace basil */
