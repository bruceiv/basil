#ifndef _DFS_HPP_
#define _DFS_HPP_

#include <vector>

#include "basilCommon.hpp"
#include "lrs/lrs.hpp"

namespace basil {
	
	/** Options for DFS algorithm.
	 *  The default constructor initializes the data members to their default 
	 *  values, while the methods can be used to modify the options.
	 */
	struct dfs_opts {
		
		/** Default constructor - sets all options to their default value */
		dfs_opts() : /*gramVec(true),*/ showAllDicts(false) {}
		
//		/** Activates (or deactivates) the gramVec option */
//		dfs_opts& usesGramVec(bool opt = true) 
//			{ gramVec = opt; return *this; }
		
		/** Activates (or deactivates) the showsAllDicts option */
		dfs_opts& showsAllDicts(bool opt = true) 
			{ showAllDicts = opt; return *this; }
		
//		/** gramVec option copied from symbal (TODO significance?) */
//		bool gramVec;
		/** show all dictionaries as they are generated */
		bool showAllDicts;
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
				: l(m), g(g), opts(opts) {}
		
		/** Perform the DFS. */
		void doDfs();
	
	private:
		
		////////////////////////////////////////////////////////////////////////
		// Typedefs for internal data types
		////////////////////////////////////////////////////////////////////////
		
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
			long index;
		};
		typedef shared_ptr<cobasis_invariants> cobasis_invariants_ptr;
		
		/** Vertex representation. */
		struct vertex_representation {
			
			vertex_representation(index_list& inc, coordinates& coords, 
					mpz_class& det) : inc(inc), coords(coords), det(det) {}
			
			index_list inc;
			coordinates coords;
			mpz_class det;
		};
		typedef shared_ptr<vertex_representation> vertex_representation_ptr;
		
		
		/** Adds a cobasis to the global list.
		 *  @param cob		The cobasis to add
		 */
		void addCobasis(cobasis_invariants_ptr cob);
		
		/** Adds a vertex to the global list
		 *  @param rep		The representation of the vertex
		 */
		void addVertex(vertex_representation_ptr rep);
		
		/** Gets the cobasis invariants for a given cobasis and coordinates
		 *  @param cob		The cobasis
		 *  @param coords	The vector coordinates
		 *  @return the invariants for the parameters
		 */
		cobasis_invariants_ptr cobasisInvariants(cobasis_ptr cob, 
				coordinates_ptr coords);
		
		/** Find the first basis for the DFS */
		void dfsFirstBasis();
		
		/** Initializes algorithm globals */
		void initGlobals();
		
		/** Gets the vertex representation for the given cobasis and 
		 *  coordinates.
		 *  @param cob		The cobasis
		 *  @param coords	The coordinates
		 *  @return the vertex representation for the parameters
		 */
		vertex_representation_ptr vertexRep(cobasis_ptr cob, 
			coordinates_ptr coords);
		
		////////////////////////////////////////////////////////////////////////
		// Initialization-time globals
		////////////////////////////////////////////////////////////////////////
		
		/** LRS wrapper for this DFS */
		lrs::lrs l;
		/** Permutation group used for this DFS */
		permutation_group& g;
		/** Options for controlling the DFS algorithm */
		dfs_opts opts;
		
		////////////////////////////////////////////////////////////////////////
		// Algorithm data
		////////////////////////////////////////////////////////////////////////
		
		/** How many bases have been found */
		long basisCount;
		/** The first cobasis found */
		cobasis_invariants_ptr initialCobasis;
		
	}; /* class dfs */

} /* namespace basil */
#endif /* _DFS_HPP_ */