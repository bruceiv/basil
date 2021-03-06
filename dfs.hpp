#ifndef _DFS_HPP_
#define _DFS_HPP_

/** Main Basil depth-first-search algorithm.
 *  NOTE: Will not include along with dfsp.hpp
 *
 *  @author Aaron Moss
 */

/*  Copyright: Aaron Moss, 2012, moss.aaron@unb.ca  */

/*  This file is part of Basil.

    Basil is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Basil is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Basil.  If not, see <http://www.gnu.org/licenses/>.  */

#include <ctime>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/functional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <gmpxx.h>

#ifdef BAS_WALLTIME
#include <sys/time.h>
#endif /* BAS_WALLTIME */

#include "basil.hpp"
#include "dfs_types.hpp"
#include "fund_domain.hpp"
#include "gram.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"


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
				: assumesNoSymmetry(false), 
				basisLimit(std::numeric_limits<unsigned long>::max()), 
				cacheSize(1000), dualFacetTrick(true), firstCobasis(),
				fundDomainLimit(0), gramVec(true), debugGram(false),
				lexOnly(false), lrs_o(), out(&std::cout), printBasis(0),
				printNew(false), printRay(0), printVertex(0),
				showsAllDicts(false), stabSearch(false) {}
		
		/** Sets (or unsets) the aRepresentation option */
		dfs_opts& inARepresentation(bool opt = true)
			{ aRepresentation = opt; return *this; }
		
		/** Activates (or deactivates) the assumesNoSymmetry option */
		dfs_opts& assumeNoSymmetry(bool opt = true)
			{ assumesNoSymmetry = opt; return *this; }
		
		/** Sets the maximum number of bases to be considered */
		dfs_opts& withBasisLimit(unsigned long lim)
			{ basisLimit = lim; return *this; }
		
		/** Sets the size of the cobasis cache */
		dfs_opts& withCacheSize(long size)
			{ cacheSize = size; return *this; }
		
		/** Deactivates (or activates) the dualFacetTrick option */
		dfs_opts& noDualFacetTrick(bool opt = true) 
			{ dualFacetTrick = !opt; return *this; }
		
		/** Sets the initial cobasis to DFS from */
		dfs_opts& withFirstCobasis(shared_ptr<lrs::index_set>& ptr) 
			{ firstCobasis = ptr; return *this; }
		
		dfs_opts& withFundDomainLimit(unsigned long lim)
			{ fundDomainLimit = lim; return *this; }

		/** Deactivates (or activates) the gramVec option */
		dfs_opts& noGramVec(bool opt = true) 
			{ gramVec = !opt; return *this; }
		
		/** Activates (or deactivates) the debugGram option */
		dfs_opts& doDebugGram(bool opt = true) 
			{ debugGram = opt; return *this; }
		
		/** Activates (or deactivates) the lexOnly option */
		dfs_opts& withLexOnly(bool opt = true)
			{ lexOnly = opt; return *this; }
		
		/** Sets the output stream for this DFS. Also sets the output stream 
		 *  for the associated LRS instance. */
		dfs_opts& withOutput(std::ostream& o) 
			{ out = &o; lrs_o.withOutput(o); return *this; }
		
		/** Gets the output stream for this DFS */
		std::ostream& output()
			{ return *out; }
		
		/** Convenience for print{Basis,Ray,Vertex} */
		dfs_opts& printAt(long n)
			{ printBasis = n; printRay = n; printVertex = n; return *this; }
		
		/** Sets the basis printing interval */
		dfs_opts& printBasisAt(long n)
			{ printBasis = n; return *this; }
		
		/** Activates (or deactivates) the printNew option */
		dfs_opts& doPrintNew(bool opt = true)
			{ printNew = opt; return *this; }
		
		/** Activates (or deactivates) the trace printing option */
		dfs_opts& doPrintTrace(bool opt = true)
			{ printTrace = opt; return *this; }
		
		/** Sets the ray printing interval */
		dfs_opts& printRayAt(long n)
			{ printRay = n; return *this; }
		
		/** Sets the vertex printing interval */
		dfs_opts& printVertexAt(long n)
			{ printVertex = n; return *this; }
		
		/** Activates (or deactivates) the showsAllDicts option */
		dfs_opts& showAllDicts(bool opt = true) 
			{ showsAllDicts = opt; return *this; }
		
		/** Activates (or deactivates) the stabSearch option */
		dfs_opts& useStabSearch(bool opt = true)
			{ stabSearch = opt; return *this; }
		
		/** Sets (or unsets) the V-representation flag  */
		dfs_opts& inVRepresentation(bool opt = true)
			{ lrs_o.inVRepresentation(opt); return *this; }
		
		
		/** Specify the input is an arrangement, and should use the arrangement 
		 *  pivoting rule */
		bool aRepresentation;
		/** assumes the given polytope is asymmetric [false]. This is primarily 
		 *  a debugging option */
		bool assumesNoSymmetry;
		/** maximum number of bases to consider
		 *  [numeric_limits\<unsigned long\>::max()] */
		unsigned long basisLimit;
		/** size of the seen cobasis lookup cache [1000] */
		long cacheSize;
		/** use the dual facet trick [true] */
		bool dualFacetTrick;
		/** cobasis to start the search from [null]. If this is not set, the 
		 *  DFS algorithm will find the first cobasis itself */
		shared_ptr<lrs::index_set> firstCobasis;
		/** maximum number of constraints to add to the fundamental domain [0]
		 */
		unsigned long fundDomainLimit;
		/** Use gram vector hashing to shrink the symmetry search space [true] 
		 */
		bool gramVec;
		/** Print gram vectors for new cobases/vertices that are printed 
		 *  [false] */
		bool debugGram;
		/** lexically based pivots only [false]. This is bad and breaks the 
		 *  algorithm, don't use it */
		bool lexOnly;
		/** options for LRS (modifier methods are cloned into here) */
		lrs::lrs_opts lrs_o;
		/** output stream for this algorithm [standard output]. */
		std::ostream* out;
		/** number of new (up to symmetry) cobases to print a progress report 
		 *  after; if 0, will never print [0]. */
		long printBasis;
		/** print the new {cobasis,ray,vertex} when seen */
		bool printNew;
		/** trace the full execution of the DFS */
		bool printTrace;
		/** number of new (up to symmetry) rays to print a progress report 
		 *  after; if 0, will never print [0]. */
		long printRay;
		/** number of new (up to symmetry) vertices to print a progress report 
		 *  after; if 0, will never print [0]. */
		long printVertex;
		/** show all dictionaries as they are generated [false] */
		bool showsAllDicts;
		/** search for cobasis symmetries using set stabilizers rather than the 
		 *  full symmetry group [false]. Not reccommended, stabilizer 
		 *  computation costs more than it saves */
		bool stabSearch;
	}; /* struct dfs_opts */
	
	/** Stateful wrapper class for DFS algorithm. */
	class dfs {
	private:
		
		////////////////////////////////////////////////////////////////////////
		// Typedefs for internal data types
		////////////////////////////////////////////////////////////////////////
		
		typedef lrs::index_set_iter index_set_iter;
		typedef lrs::index_set_hash index_set_hash;
		
		typedef lrs::vector_mpz vector_mpz;
		typedef shared_ptr<vector_mpz> vector_mpz_ptr;
		typedef lrs::vector_mpq_hash coordinates_hash;
		typedef lrs::matrix_mpq_hash matrix_hash;
		
		typedef lrs::cobasis cobasis;
		typedef shared_ptr<cobasis> cobasis_ptr;
		
		typedef std::pair<ind, ind> index_pair;
	
	public:
		
		/** Set up a DFS on the given matrix, with the given permuation group.
		 *  @param m		The matrix to DFS on
		 *  @param lin		The set of indices that are linearities
		 *  @param g		The permutation group of the matrix
		 *  @param gram		The gram matrix for the constraints (may be empty 
		 *  				if o.gramVec == false)
		 *  @param opts		The options for this DFS (default values if not 
		 * 					provided)
		 */
		dfs(matrix& m, index_set& lin, permutation_group& g, gram_matrix& gram, 
				dfs_opts o = dfs_opts());
		
		
		/** Perform the DFS.
		 *  @return did the DFS run to completion, or halt because of too many 
		 *  	bases?
		 */
		bool doDfs();
		
		////////////////////////////////////////////////////////////////////////
		// Query methods for after completion of doDfs()
		////////////////////////////////////////////////////////////////////////
		
		/** @return representatives of each of the orbits of the cobases */
		cobasis_map const& getBasisOrbits() const;
		
		/** @return the sum of the degrees of the basis orbit representatives */
		uind getTotalBasisDegree() const;

		/** @return the dimension of the polytope */
		ind getDimension() const;
		
		/** @return the initial cobasis for the DFS */
		index_set getInitialCobasis() const;
		
		/** @return did the DFS complete (true) or terminate due to too many 
		 *  	bases (false) */
		bool isFinished() const;
		
		/** @return the generated fundamental domain */
		fund_domain const& getFundamentalDomain() const;

		/** @return representatives of each of the orbits of the extreme rays */
		coordinates_map const& getRayOrbits() const;
		
		/** @return the running time of the algorithm, in milliseconds */
		std::clock_t getRunningTime() const;
		
#ifdef BAS_WALLTIME
		/** @return the wall time the algorithm took, in milliseconds */
		long getWallTime() const;
#endif /* BAS_WALLTIME */

		/** @return the symmetry group used in the DFS */
		permutation_group const& getSymmetryGroup() const;
		
		/** @return representatives of each of the orbits of the vertices */
		coordinates_map const& getVertexOrbits() const;
		
		/** @return the gram matrix. (This method is of primary use for 
		 *  debugging the gram vectors.) */
		gram_matrix const& getGramMat() const;
		
		/** @return the map of gram vectors to cobases. (This method is of 
		 *  primary use for debugging the gram vectors.) */
		cobasis_gram_map const& getCobasisGramMap() const;
		
		/** @return the map of gram vectors to vertices. (This method is of 
		 *  primary use for debugging the gram vectors.) */
		vertex_gram_map const& getVertexGramMap() const;
		
	private:
		
		////////////////////////////////////////////////////////////////////////
		// More internal data structures
		////////////////////////////////////////////////////////////////////////
		
		/** Representation of a pivot */
		struct pivot {
			
			pivot(index_set cob, ind leave, ind enter) 
					: cob(cob), leave(leave), enter(enter) {}
			
			/** cobasis before pivot */
			index_set cob;
			/** leaving index */
			ind leave;
			/** entering index */
			ind enter;
		};
		typedef shared_ptr<pivot> pivot_ptr;
		
		typedef 
			std::pair<cobasis_gram_map::const_iterator, 
				cobasis_gram_map::const_iterator> 
			cobasis_gram_range;
		
		typedef
			std::pair<vertex_gram_map::const_iterator, 
				vertex_gram_map::const_iterator>
			vertex_gram_range;
		
		/** Transformation of index_set_iter to provide input to PermLib 
		 *  properly. The public interfaces to PermLib are zero-indexed, 
		 *  whereas the I/O operators, rather inconsistently, are one-indexed; 
		 *  this type provides a suitable wrapper for conversion of a 
		 *  one-indexed index_set_iter into a zero-indexed iterator to provide 
		 *  input to PermLib.
		 */
		typedef 
			boost::transform_iterator< 
				boost::binder2nd<
					std::minus<index_set_iter::value_type>
				>,
				index_set_iter
			> 
			pl_index_set_iter;
		
		/** Transform of lrs::begin(s) iterator to PermLib zero-index. */
		pl_index_set_iter plBegin(index_set& s) {
			return boost::make_transform_iterator(
					lrs::begin(s), 
					boost::bind2nd(std::minus<index_set_iter::value_type>(), 1)
			);
		}
		
		/** Transform of lrs::end(s) iterator to PermLib zero-index. */
		pl_index_set_iter plEnd(index_set& s) {
			return boost::make_transform_iterator(
					lrs::end(s), 
					boost::bind2nd(std::minus<index_set_iter::value_type>(), 1)
			);
		}
		
		/** Adds a cobasis to the global map of cobasis representatives. Note 
		 *  that the vertex data pointer given should already be present in the 
		 *  global map of vertex representatives (added by another cobasis), or 
		 *  addVertex() should be used (which calls this method).
		 *  @param cob		The cobasis to map to the vertex
		 *  @param dat		The data for this vertex
		 */
		void addCobasis(index_set const& cob, vertex_data_ptr dat);
		
		/** Adds a vertex to the global map of vertex representatives.
		 *  This will add all cobases defined for this vertex to the global map 
		 *  of cobasis representatives as well.
		 *  @param dat		The data for this vertex
		 */
		void addVertex(vertex_data_ptr dat);
		
		/** Find the first vertex for the DFS */
		index_set dfsFirstBasis();
		
		/** DFS from a starting cobasis. Caller is responsible for pivoting to 
		 *  the correct dictionary before calling. Will set the hitMaxBasis 
		 *  variable appropriately before returning.
		 *  @param root		The root cobasis to DFS from
		 *  @return did the DFS run to completion, or halt because of too many 
		 *  	bases? (will set hitMaxBasis to the negation of this before 
		 *  	returning)
		 */
		bool dfsFromRoot(index_set& root);
		
		/** Gets the gram vector for an incidence set.
		 *  @param inc		The incidence set to take the gram vector for
		 *  @return the gram matrix, restricted to this incidence set in row 
		 *  		and column indices, then sorted
		 */
		gram_matrix fastGramVec(index_set inc);
		
		/** Looks for symmetries between a given cobasis and a list of 
		 *  candidate cobases.
		 *  @return true if a symmetry is found, false otherwise
		 */
		bool findSymmetry(index_set find, index_set_list list);
		
		/** Finds the rays in the current dictionary. */
		void getRays();
		
		/** Checks if a cobasis has been seen before.
		 *  @param cob		The cobasis to check
		 *  @param dat		The invariant data for this cobasis
		 *  @return true for likely seen, false for likely not
		 */
		bool isNewCobasis(index_set cob, vertex_data_ptr dat);
		
		/** Gets the canonical ray for each ray in a known orbit.
		 *  @param rep		The ray to get the orbit representative of
		 *  @return a pointer to the ray representative of this ray's orbit, or 
		 * 		a null pointer if there is none such.
		 */
		vertex_data_ptr knownRay(vertex_data_ptr rep);
		
		/** Gets the canonical vertex for each vertex in a known orbit.
		 *  @param rep		The vertex to get the orbit representative of
		 *  @return a pointer to the vertex representative of this vertex's 
		 * 		orbit, or a null pointer if there is none such.
		 */
		vertex_data_ptr knownVertex(vertex_data_ptr rep);
		
		/** Gets the cobases whose invariants match the given one.
		 *  @param cob		The cobasis to match
		 *  @param dat		The invariant data to match
		 *  @return a list of cobases with matching invariants.
		 */
		index_set_list matchingCobasisInvariants(index_set cob, 
												 vertex_data_ptr dat);
		
		/** Gets the vertices whose invariants match the given one.
		 *  @param rep		The vertex data to match
		 *  @return a list of vertices with matching invariants.
		 */
		vertex_data_list matchingInvariants(vertex_data_ptr rep);
		
		/** Check the invariants of two vertex data for cobasis symmetry
		 *  @param a			The first vertex to check
		 *  @param b			The second vertex to check
		 *  @return true if the invariants of the two vertex data match. 
		 */
		bool cobasisInvariantsMatch(vertex_data const& a, vertex_data const& b);
		
		/** Check the invariants of two vertex data for vertex symmetry
		 *  @param a			The first vertex to check
		 *  @param b			The second vertex to check
		 *  @return true if the invariants of the two vertex data match. 
		 */
		bool invariantsMatch(vertex_data const& a, vertex_data const& b);
		
		/** Add new edges to the search stack.
		 *  @param oldCob	The cobasis to search for adjacent edges
		 */
		void pushNewEdges(index_set& oldCob);
		
		/** Combines the cobasis and coordinate data into a vertex_data object
		 *  @param cob		The cobasis data
		 *  @param coords	The ray coordinates
		 *  @return the vertex data from the parameters
		 */
		vertex_data_ptr rayData(cobasis_ptr cob, vector_mpz_ptr coords);
		
		/** Combines the cobasis and coordinate data into a vertex_data object
		 *  @param cob		The cobasis data
		 *  @param coords	The vertex coordinates (un-rationalized)
		 *  @return the vertex data from the parameters
		 */
		vertex_data_ptr vertexData(cobasis_ptr cob, vector_mpz_ptr coords);
		
		////////////////////////////////////////////////////////////////////////
		// Initialization-time globals
		////////////////////////////////////////////////////////////////////////
		
		/** LRS wrapper for this DFS */
		lrs::lrs l;
		/** Permutation group used for this DFS */
		permutation_group& g;
		/** Original constraint matrix */
		matrix& m;
		/** Options for controlling the DFS algorithm */
		dfs_opts opts;
		/** Dimension of the problem */
		ind dim;
		/** number of rows in the problem */
		ind rows;
		/** matrix of pre-computed inner product representatives, for gram 
		 *  vectors */
		gram_matrix gramMat;
		
		////////////////////////////////////////////////////////////////////////
		// Algorithm data
		////////////////////////////////////////////////////////////////////////
		
		/** A constant list of all the indices */
		index_set allIndices;
		/** Cache of recently seen cobases */
		lru::cache<index_set, index_set_hash> cobasisCache;
		/** Lookup cobases by gram vector */
		cobasis_gram_map cobasisGramMap;
		/** Global map of seen cobases, up to symmetry */
		cobasis_map basisOrbits;
		/** Search queue for cobases */
		std::deque<index_set> cobasisQueue;
		/** Sum of the degrees of the cobasis orbit representatives */
		uind totalBasisDegree;
		/** Temporary to store time diffs into. Will hold total running time on 
		 *  algorithm completion. */
		std::clock_t diff_time;
#ifdef BAS_WALLTIME
		/** Temporary to store time diffs into. Will hold total running time on
		 *  algorithm completion. Uses wall time instead of CPU time. */
		struct timeval wall_end_time;
#endif /* BAS_WALLTIME */
		/** If the basis count hit the maximum count */
		bool hitMaxBasis;
		/** The first cobasis found */
		index_set initialCobasis;
		/** The fundamental domain of the problem, iteratively constructed */
		fund_domain fundDomain;
		/** Backtracking stack. */
		std::deque<index_pair> pathStack;
		/** representatives of each orbit (of rays) */
		coordinates_map rayOrbits;
		/** The true dimension of the polytope */
		ind realDim;
		/** the time at which the algorithm was started */
		std::clock_t start_time;
#ifdef BAS_WALLTIME
		/* The wall time at which the algorithm was started */
		struct timeval wall_start_time;
#endif /* BAS_WALLTIME */
		/** representatives of each orbit (of vertices) */
		coordinates_map vertexOrbits;
		/** Lookup vertices by gram vector */
		vertex_gram_map vertexGramMap;
		/** Pivots in the working stack */
		std::deque<pivot> workStack;
		
		////////////////////////////////////////////////////////////////////////
		// Time-related functions and values
		////////////////////////////////////////////////////////////////////////
		
		/** Number of clock ticks for millisecond */
		static std::clock_t const clocks_per_ms = CLOCKS_PER_SEC / 1000;
		
		/** Update and return the current running time. */
		std::clock_t currentTime() {
			diff_time = std::clock() - start_time;
			return diff_time / clocks_per_ms;
		}
		
	}; /* class dfs */

} /* namespace basil */
#endif /* _DFS_HPP_ */
