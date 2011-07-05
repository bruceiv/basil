#ifndef _DFS_HPP_
#define _DFS_HPP_

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

#include <permlib/common.h> //because the PermLib author didn't ...
#include <permlib/bsgs.h>
#include <permlib/permutation.h>
#include <permlib/transversal/schreier_tree_transversal.h>

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"


namespace basil {
	
	////////////////////////////////////////////////////////////////////////////
	//
	//  Imports and typedefs for use in Basil
	//
	////////////////////////////////////////////////////////////////////////////
	
	/** import STL string into this namespace */
	using std::string;
	
	/** import boost shared pointer into this namespace */
	using boost::shared_ptr;
	

	/** matrix type */
	typedef 
		lrs::matrix
		matrix;
	typedef
		shared_ptr<matrix>
		matrix_ptr;
	
	/** typesafe index into matrix */
	typedef 
		lrs::ind
		ind;
	/** unsigned version of ind */
	typedef
		lrs::uind
		uind;
	
	/** permutation type */
	typedef 
		permlib::Permutation 
		permutation;
	typedef
		shared_ptr<permutation>
		permutation_ptr;
		
	/** permutation tree traversal type */
	typedef 
		permlib::SchreierTreeTransversal<permutation>
		permutation_transversal;
	typedef
		shared_ptr<permutation_transversal>
		permutation_transversal_ptr;
	
	/** permutation group type */
	typedef 
		permlib::BSGS<permutation, permutation_transversal> 
		permutation_group;
	typedef
		shared_ptr<permutation_group>
		permutation_group_ptr;
	
	/** list of permutation type */
	typedef
		typename permutation_group::PERMlist
		permutation_list;
	
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
				gramVec(true), lexOnly(false), lrs_o(), out(&std::cout), 
				printBasis(0), printNew(false), printRay(0), printVertex(0), 
				showsAllDicts(false) {}
		
		
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
		
		/** Deactivates (or activates) the gramVec option */
		dfs_opts& noGramVec(bool opt = true) 
			{ gramVec = !opt; return *this; }
		
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
		
		/** Sets the ray printing interval */
		dfs_opts& printRayAt(long n)
			{ printRay = n; return *this; }
		
		/** Sets the vertex printing interval */
		dfs_opts& printVertexAt(long n)
			{ printVertex = n; return *this; }
		
		/** Activates (or deactivates) the showsAllDicts option */
		dfs_opts& showAllDicts(bool opt = true) 
			{ showsAllDicts = opt; return *this; }
		
		/** Sets (or unsets) the V-representation flag  */
		dfs_opts& inVRepresentation(bool opt = true)
			{ lrs_o.inVRepresentation(opt); return *this; }
		
		
		/** assumes the given polytope is asymmetric [false]. This is primarily 
		 *  a debugging option */
		bool assumesNoSymmetry;
		/** maximum number of bases to consider [numeric_limits\<long\>::max()] 
		 */
		unsigned long basisLimit;
		/** size of the seen cobasis lookup cache [1000] */
		long cacheSize;
		/** use the dual facet trick [true] */
		bool dualFacetTrick;
		/** cobasis to start the search from [null]. If this is not set, the 
		 *  DFS algorithm will find the first cobasis itself */
		shared_ptr<lrs::index_set> firstCobasis;
		/** Use gram vector hashing to shrink the cobasis search space [true] */
		bool gramVec;
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
		/** number of new (up to symmetry) rays to print a progress report 
		 *  after; if 0, will never print [0]. */
		long printRay;
		/** number of new (up to symmetry) vertices to print a progress report 
		 *  after; if 0, will never print [0]. */
		long printVertex;
		/** show all dictionaries as they are generated [false] */
		bool showsAllDicts;
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
		
		typedef lrs::cobasis cobasis;
		typedef shared_ptr<cobasis> cobasis_ptr;
		
		typedef std::pair<ind, ind> index_pair;
	
	public:
		
		////////////////////////////////////////////////////////////////////////
		// Typedefs for external data types
		////////////////////////////////////////////////////////////////////////
		
		typedef lrs::vector_mpq coordinates;
		
		typedef lrs::index_set index_set;
		typedef shared_ptr<index_set> index_set_ptr;
		typedef std::vector<index_set> index_set_list;
		
		/** Joint vertex-cobasis storage */
		struct vertex_data {
			
			/** Single-cobasis constructor. Initializes all fields as you would 
			 *  think, where cobs is set up to be a set initially including 
			 *  only cob. 
			 */
			vertex_data(coordinates coords, index_set inc, index_set cob, 
					mpz_class det) : coords(coords), inc(inc), cobs(), 
					det(det) {
				cobs.insert(cob);
			}
			
			/** Multiple-cobasis constructor. Initializes all fields to the 
			 *  given values */
			vertex_data(coordinates coords, index_set inc, 
					std::set<index_set> cobs, mpz_class det) : coords(coords), 
					inc(inc), cobs(cobs), det(det) { }
			
			/* Key data */
			/** Coordinates of the vertex */
			coordinates coords;
			/** Set of incident cobasis indices */
			index_set inc;
			/** Set of cobases for this vertex */
			std::set<index_set> cobs;
			
			/* Invariants */
			/** determinant TODO (?) */
			mpz_class det;
		};
		typedef shared_ptr<vertex_data> vertex_data_ptr;
		
		/** map of vertex coordinates to a vertex data pointer */
		typedef 
			boost::unordered_map<coordinates, vertex_data_ptr, coordinates_hash>
			coordinates_map;
		/** map of a cobasis to a vertex data pointer */
		typedef
			boost::unordered_map<index_set, vertex_data_ptr, index_set_hash>
			cobasis_map;
		
		
		/** Set up a DFS on the given matrix, with the given permuation group.
		 *  @param m		The matrix to DFS on
		 *  @param g		The permutation group of the matrix
		 *  @param opts		The options for this DFS (default values if not 
		 * 					provided)
		 */
		dfs(matrix& m, permutation_group& g, dfs_opts o = dfs_opts()) 
				: l(m, o.lrs_o), g(g), opts(o) { 
			dim = m.d();
			rows = m.n();
			
			/* set up algorithm globals */
			initGlobals();
		}
		
		
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
		
		/** @return the dimension of the polytope */
		ind getDimension() const;
		
		/** @return the initial cobasis for the DFS */
		index_set getInitialCobasis() const;
		
		/** @return did the DFS complete (true) or terminate due to too many 
		 *  	bases (false) */
		bool isFinished() const;
		
		/** @return representatives of each of the orbits of the extreme rays */
		coordinates_map const& getRayOrbits() const;
		
		/** @return the running time of the algorithm, in milliseconds */
		std::clock_t getRunningTime() const;
		
		/** @return the symmetry group used in the DFS */
		permutation_group const& getSymmetryGroup() const;
		
		/** @return representatives of each of the orbits of the vertices */
		coordinates_map const& getVertexOrbits() const;
		
		
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
		
		/** Finds the rays in the current dictionary. */
		void getRays();
		
		/** Looks for symmetries between a given cobasis and a list of 
		 *  candidate cobases.
		 *  @return true if a symmetry is found, false otherwise
		 */
		bool findSymmetry(index_set find, index_set_list list);
		
		/** Initializes algorithm globals */
		void initGlobals();
		
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
		 *  @param dat		The invariant data to match
		 *  @return a list of cobases with matching invariants.
		 */
		index_set_list matchingInvariants(vertex_data_ptr dat);
		
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
		index_set allIndices;
		/** Cache of recently seen cobases */
		lru::cache<index_set, index_set_hash> cobasisCache;
		/** Global map of seen cobases, up to symmetry */
		cobasis_map basisOrbits;
		/** Search queue for cobases */
		std::deque<index_set> cobasisQueue;
		/** Temporary to store time diffs into. Will hold total running time on 
		 *  algorithm completion. */
		std::clock_t diff_time;
		/** If the basis count hit the maximum count */
		bool hitMaxBasis;
		/** The first cobasis found */
		index_set initialCobasis;
		/** Backtracking stack. */
		std::deque<index_pair> pathStack;
		/** representatives of each orbit (of rays) */
		coordinates_map rayOrbits;
		/** The true dimension of the polytope */
		ind realDim;
		/** the time at which the algorithm was started */
		std::clock_t start_time;
		/** representatives of each orbit (of vertices) */
		coordinates_map vertexOrbits;
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