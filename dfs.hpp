#ifndef _DFS_HPP_
#define _DFS_HPP_

#include <deque>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "basilCommon.hpp"
#include "lruCache.hpp"

#include "lrs/lrs.hpp"


namespace basil {
	
	/** Exception thrown for unexpected circumstances in the DFS algorithm.
	 *  The what string will describe the error.
	 */
	class dfs_error : public std::runtime_error {
	public:
		dfs_error(std::string const& whatArg) : runtime_error(whatArg) {}
	};
	
	/** Options for DFS algorithm.
	 *  The default constructor initializes the data members to their default 
	 *  values, while the methods can be used to modify the options.
	 */
	struct dfs_opts {
		
		/** Default constructor - sets all options to their default value */
		dfs_opts() 
				: cacheSize(1000), dualFacetTrick(true), 
				showsAllDicts(false) {}
		
		/** Sets the size of the cobasis cache */
		dfs_opts& withCacheSize(long size)
			{ cacheSize = size; return *this; }
		
		/** Deactivates (or activates) the dualFacetTrick option */
		dfs_opts& noDualFacetTrick(bool opt = true) 
			{ dualFacetTrick = !opt; return *this; }
		
		/** Activates (or deactivates) the showsAllDicts option */
		dfs_opts& showAllDicts(bool opt = true) 
			{ showsAllDicts = opt; return *this; }
		
		/** size of the seen cobasis lookup cache [1000] */
		long cacheSize;
		/** use the dual facet trick [true] */
		bool dualFacetTrick;
		/** show all dictionaries as they are generated [false] */
		bool showsAllDicts;
	}; /* struct dfs_opts */
	
	/** Stateful wrapper class for DFS algorithm. */
	class dfs {
	public:
		/** Set up a DFS on the given matrix, with the given permuation group.
		 *  @param m		The matrix to DFS on
		 *  @param g		The permutation group of the matrix
		 *  @param opts		The options for this DFS (default values if not 
		 * 					provided)
		 */
		dfs(matrix& m, permutation_group& g, dfs_opts opts = dfs_opts()) 
				: l(m), g(g), opts(opts), dim(m.d), rows(m.n) {}
		
		/** Perform the DFS. */
		void doDfs();
	
	private:
		
		////////////////////////////////////////////////////////////////////////
		// Typedefs for internal data types
		////////////////////////////////////////////////////////////////////////
		
		typedef lrs::ind ind;
		
		typedef lrs::index_list index_list;
		typedef shared_ptr<index_list> index_list_ptr;
		
		typedef lrs::cobasis cobasis;
		typedef shared_ptr<cobasis> cobasis_ptr;
		
		typedef lrs::vector_mpz coordinates;
		typedef shared_ptr<coordinates> coordinates_ptr;
		
		/** Invariants of a cobasis */
		struct cobasis_invariants {
			
			cobasis_invariants(index_list& cob, index_list& extraInc, 
					coordinates& coords, mpz_class& det)
					: cob(cob), extraInc(extraInc), coords(coords), det(det)
					{}
			
			index_list cob;
			index_list extraInc;
			coordinates coords;
			mpz_class det;
			ind index;
		};
		typedef shared_ptr<cobasis_invariants> cobasis_invariants_ptr;
		typedef std::vector<cobasis_invariants_ptr> cobasis_invariants_list;
		
		/** Vertex representation. */
		struct vertex_rep {
			
			vertex_rep(index_list& inc, coordinates& coords, 
					mpz_class& det) : inc(inc), coords(coords), det(det) {}
			
			/* NOTE: should be kept sorted */
			index_list inc;
			coordinates coords;
			mpz_class det;
		};
		typedef shared_ptr<vertex_rep> vertex_rep_ptr;
		
		/** Representation of a pivot */
		struct pivot {
			
			pivot(index_list& cob, ind leave, ind enter) 
					: cob(cob), leave(leave), enter(enter) {}
			
			/** cobasis before pivot */
			index_list cob;
			/** leaving index */
			ind leave;
			/** entering index */
			ind enter;
		};
		typedef shared_ptr<pivot> pivot_ptr;
		
		
		/** Adds a cobasis to the global list.
		 *  @param cob		The cobasis to add
		 */
		void addCobasis(cobasis_invariants_ptr cob);
		
		/** Adds a vertex to the global list
		 *  @param rep		The representation of the vertex
		 */
		void addVertex(vertex_rep_ptr rep);
		
		/** Gets the cobasis invariants for a given cobasis and coordinates
		 *  @param cob		The cobasis
		 *  @param coords	The vector coordinates
		 *  @return the invariants for the parameters
		 */
		cobasis_invariants_ptr cobasisInvariants(cobasis_ptr cob, 
				coordinates_ptr coords);
		
		/** Find the first basis for the DFS */
		cobasis_ptr dfsFirstBasis();
		
		/** DFS from a starting cobasis. Caller is responsible for pivoting to 
		 *  the correct dictionary before calling.
		 *  @param root		The root cobasis to DFS from
		 */
		void dfsFromRoot(index_list& root);
		
		/** Finds the rays in the current dictionary. */
		void getRays();
		
		/** Looks for symmetries between a given cobasis and a list of 
		 *  candidate cobases.
		 *  @return true if a symmetry is found, false otherwise
		 */
		bool findSymmetry(cobasis_invariants_ptr rep, 
						  cobasis_invariants_list list);
		
		/** Initializes algorithm globals */
		void initGlobals();
		
		/** Checks if a given cobasis has been seen before. 
		 *  @param rep		The cobasis to check
		 *  @return true for likely seen, false for likely not
		 */
		bool isNewCobasis(cobasis_invariants_ptr rep);
		
		/** Gets the canonical ray for each ray in a known orbit.
		 *  @param rep		The ray to get the orbit representative of
		 *  @return a pointer to the ray representative of this ray's orbit, or 
		 * 		a null pointer if there is none such.
		 */
		vertex_rep_ptr knownRay(vertex_rep_ptr rep);
		
		/** Gets the canonical vertex for each vertex in a known orbit.
		 *  @param rep		The vertex to get the orbit representative of
		 *  @return a pointer to the vertex representative of this vertex's 
		 * 		orbit, or a null pointer if there is none such.
		 */
		vertex_rep_ptr knownVertex(vertex_rep_ptr rep);
		
		/** Gets the cobases whose invariants match the given one.
		 *  @return a list of cobases with matching invariants.
		 */
		cobasis_invariants_list matchingInvariants(cobasis_invariants_ptr rep);
		
		/** Add new edges to the search stack.
		 *  @param oldCob	The cobasis to search for adjacent edges
		 */
		void pushNewEdges(index_list& oldCob);
		
		/** Gets the ray representation for the given cobasis and 
		 *  coordinates.
		 *  @param cob		The cobasis
		 *  @param coords	The coordinates
		 *  @return the ray representation for the parameters
		 */
		vertex_rep_ptr rayRep(cobasis_ptr cob, coordinates_ptr coords);
		
		/** Gets the vertex representation for the given cobasis and 
		 *  coordinates.
		 *  @param cob		The cobasis
		 *  @param coords	The coordinates
		 *  @return the vertex representation for the parameters
		 */
		vertex_rep_ptr vertexRep(cobasis_ptr cob, coordinates_ptr coords);
		
		////////////////////////////////////////////////////////////////////////
		// Initialization-time globals
		////////////////////////////////////////////////////////////////////////
		
		/** Dimension of the problem */
		ind dim;
		/** LRS wrapper for this DFS */
		lrs::lrs l;
		/** Permutation group used for this DFS */
		permutation_group& g;
		/** Options for controlling the DFS algorithm */
		dfs_opts opts;
		/** number of rows in the problem */
		ind rows;
		
		////////////////////////////////////////////////////////////////////////
		// Algorithm data
		////////////////////////////////////////////////////////////////////////
		
		/** A constant list of all the indices */
		index_list allIndices;
		/** How many bases have been found */
		ind basisCount;
		/** Cache of recently seen cobases */
		lru_cache<index_list> cobasisCache;
		/** Global list of seen cobases */
		cobasis_invariants_list cobasisList;
		/** Search queue for cobases */
		std::deque<index_list> cobasisQueue;
		/** The first cobasis found */
		cobasis_invariants_ptr initialCobasis;
		/** representatives of each orbit (of rays) */
		std::vector<vertex_rep_ptr> rayOrbits;
		/** The true dimension of the polytope */
		ind realDim;
		/** representatives of each orbit (of vertices) */
		std::vector<vertex_rep_ptr> vertexOrbits;
		/** coordinates found */
		std::set<coordinates> vertexSet;
		/** Pivots in the working stack */
		std::deque<pivot> workStack;
		
	}; /* class dfs */

} /* namespace basil */
#endif /* _DFS_HPP_ */