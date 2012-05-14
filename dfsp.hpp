#ifndef _DFS_HPP_
#define _DFS_HPP_

/* NOTE: Do not include with dfs.hpp - name collisions */

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

#include "basil.hpp"
#include "gram.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

#include "lru/cache.hpp"

namespace basil {
	
	/** Exception thrown for unexpected circumstances in the DFS 
	 *  algorithm. The what string will describe the error.
	 */
	class dfs_error : public std::runtime_error {
	public:
		dfs_error(std::string const& what) : runtime_error(what) {}
	};
	
	/** Options for DFS algorithm.
	 *  The default constructor initializes the data members to their 
	 *  default values, while the methods can be used to modify the 
	 *  options.
	 */
	struct dfs_opts {
		
		/** Default constructor - sets all options to their default 
		 *  value. */
		dfs_opts() 
				: assumesNoSymmetry(false), 
				basisLimit(std::numeric_limits<unsigned long>::max()), 
				cacheSize(1000), dualFacetTrick(true), firstCobasis(), 
				gramVec(true), debugGram(false), lexOnly(false), 
				lrs_o(), out(&std::cout), printBasis(0), 
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
		
		/** Deactivates (or activates) the gramVec option */
		dfs_opts& noGramVec(bool opt = true) 
			{ gramVec = !opt; return *this; }
		
		/** Activates (or deactivates) the debugGram option */
		dfs_opts& doDebugGram(bool opt = true) 
			{ debugGram = opt; return *this; }
		
		/** Activates (or deactivates) the lexOnly option */
		dfs_opts& withLexOnly(bool opt = true)
			{ lexOnly = opt; return *this; }
		
		/** Sets the output stream for this DFS. Also sets the output 
		 *  stream for the associated LRS instance. */
		dfs_opts& withOutput(std::ostream& o) 
			{ out = &o; lrs_o.withOutput(o); return *this; }
		
		/** Gets the output stream for this DFS */
		std::ostream& output()
			{ return *out; }
		
		/** Convenience for print{Basis,Ray,Vertex} */
		dfs_opts& printAt(long n) {
			printBasis = n; printRay = n; printVertex = n; 
			return *this;
		}
		
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
		
		
		/** Specify the input is an arrangement, and should use the 
		 *  arrangement pivoting rule */
		bool aRepresentation;
		/** assumes the given polytope is asymmetric [false]. This is 
		 *  primarily a debugging option */
		bool assumesNoSymmetry;
		/** maximum number of bases to consider 
		 *  [numeric_limits\<long\>::max()] */
		unsigned long basisLimit;
		/** size of the seen cobasis lookup cache [1000] */
		long cacheSize;
		/** use the dual facet trick [true] */
		bool dualFacetTrick;
		/** cobasis to start the search from [null]. If this is not 
		 *  set, the DFS algorithm will find the first cobasis itself */
		shared_ptr<lrs::index_set> firstCobasis;
		/** Use gram vector hashing to shrink the symmetry search space 
		 *  [true] */
		bool gramVec;
		/** Print gram vectors for new cobases/vertices that are 
		 *  printed [false] */
		bool debugGram;
		/** lexically based pivots only [false]. This is bad and breaks 
		 *  the algorithm, don't use it */
		bool lexOnly;
		/** options for LRS (modifier methods are cloned into here) */
		lrs::lrs_opts lrs_o;
		/** output stream for this algorithm [standard output]. */
		std::ostream* out;
		/** number of new (up to symmetry) cobases to print a progress 
		 *  report after; if 0, will never print [0]. */
		long printBasis;
		/** print the new {cobasis,ray,vertex} when seen */
		bool printNew;
		/** trace the full execution of the DFS */
		bool printTrace;
		/** number of new (up to symmetry) rays to print a progress 
		 *  report after; if 0, will never print [0]. */
		long printRay;
		/** number of new (up to symmetry) vertices to print a progress 
		 *  report after; if 0, will never print [0]. */
		long printVertex;
		/** show all dictionaries as they are generated [false] */
		bool showsAllDicts;
		/** search for cobasis symmetries using set stabilizers rather 
		 *  than the full symmetry group [false]. Not reccommended, 
		 *  stabilizer computation costs more than it saves */
		bool stabSearch;
	}; /* struct dfs_opts */
	
	/** Stateful wrapper class for DFS algorithm. */
	class dfs {
	private:

		////////////////////////////////////////////////////////////////
		// Typedefs for internal data types
		////////////////////////////////////////////////////////////////

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
		
		////////////////////////////////////////////////////////////////
		// Typedefs for external data types
		////////////////////////////////////////////////////////////////
		
		typedef lrs::vector_mpq coordinates;
		
		typedef std::vector<index_set> index_set_list;
		
		/** Joint vertex-cobasis storage */
		struct vertex_data {
			
			/** Single-cobasis constructor. Initializes all fields as 
			 *  you would think, where cobs is set up to be a set 
			 *  initially including only cob. 
			 */
			vertex_data(coordinates coords, index_set inc, 
					index_set cob, mpz_class det, gram_matrix gram) 
					: coords(coords), inc(inc), cobs(), det(det), 
					gram(gram) {
				cobs.insert(cob);
			}
			
			/** Multiple-cobasis constructor. Initializes all fields to 
			 *  the given values */
			vertex_data(coordinates coords, index_set inc, 
					std::set<index_set> cobs, mpz_class det, 
					gram_matrix gram) 
					: coords(coords), inc(inc), cobs(cobs), det(det), 
					gram(gram) { }
			
			/* Key data */
			/** Coordinates of the vertex */
			coordinates coords;
			/** Set of incident cobasis indices */
			index_set inc;
			/** Set of cobases for this vertex */
			std::set<index_set> cobs;
			
			/* Invariants */
			/** determinant */
			mpz_class det;
			/** gram matrix */
			gram_matrix gram;
		};
		typedef shared_ptr<vertex_data> vertex_data_ptr;
		typedef std::vector<vertex_data_ptr> vertex_data_list;
		
		/** map of vertex coordinates to a vertex data pointer */
		typedef 
			boost::unordered_map<
				coordinates, vertex_data_ptr, coordinates_hash>
			coordinates_map;
		/** map of a cobasis to a vertex data pointer */
		typedef
			boost::unordered_map<
				index_set, vertex_data_ptr, index_set_hash>
			cobasis_map;
		/** map of a gram vector to its cobases, and the associated vertex 
		 *  data */
		typedef
			boost::unordered_multimap<
				gram_matrix, 
				std::pair<index_set, vertex_data_ptr>, 
				gram_matrix_hash>
			cobasis_gram_map;
		/** map of a gram vector to its vertices */
		typedef
			boost::unordered_multimap<
				gram_matrix, vertex_data_ptr, gram_matrix_hash>
			vertex_gram_map;
	
	private:
		
		////////////////////////////////////////////////////////////////
		// More typedefs for internal data types
		////////////////////////////////////////////////////////////////
		
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
		
		/** A pointer to a vertex data struct, with an attatched 
		 *  boolean to indicate if the struct was previously in the 
		 *  searched structure. */
		typedef std::pair<vertex_data_ptr, bool> vertex_data_known;
		
		/** List of cobases with associated vertex data */
		typedef 
			std::vector< std::pair<index_set, vertex_data_ptr> >
			cobasis_list;
		/** List of coordinates with associated vertex data */
		typedef 
			std::vector< std::pair<coordinates, vertex_data_ptr> >
			coordinates_list;
		
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
		static pl_index_set_iter plBegin(index_set& s) {
			return boost::make_transform_iterator(
					lrs::begin(s),
					boost::bind2nd(std::minus<index_set_iter::value_type>(), 1)
			);
		}

		/** Transform of lrs::end(s) iterator to PermLib zero-index. */
		static pl_index_set_iter plEnd(index_set& s) {
			return boost::make_transform_iterator(
					lrs::end(s),
					boost::bind2nd(std::minus<index_set_iter::value_type>(), 1)
			);
		}

		////////////////////////////////////////////////////////////////
		// Thread-local data structures
		////////////////////////////////////////////////////////////////
		
		/** Thread-local information for the DFS algorithm. */
		struct explorer {
			
			/** Set up a DFS explorer on the given matrix, with the 
			 *  given permuation group.
			 *  @param m		The matrix to DFS on
			 *  @param lin		The set of indices that are linearities
			 *  @param g		The permutation group of the matrix
			 *  @param gram		The gram matrix for the constraints 
			 *  				(may be empty if o.gramVec == false)
			 */
			explorer(matrix m, index_set lin, permutation_group g,
					gram_matrix gram, dfs_opts o = dfs_opts());
			
			/** Gets the canonical ray for each ray in a known orbit.
			 *  @param rays		The set of ray orbit representatives
			 *  @param rep		The ray to get the orbit representative 
			 *  				of
			 *  @return a pointer to the ray representative of this 
			 *  		ray's orbit, or a null pointer if there is none 
			 *  		such.
			 */
			vertex_data_ptr knownRay(coordinates_map rays, 
					vertex_data_ptr rep);
			
			////////////////////////////////////////////////////////////
			// Thread-local copies of initilization time globals
			////////////////////////////////////////////////////////////
			
			/** LRS wrapper for this DFS thread */
			lrs::lrs l;
			/** Permutation group used for this DFS */
			permutation_group& g;
			/** matrix of pre-computed inner product representatives, 
			 *  for gram vectors */
			gram_matrix gramMat;
			/** Options for controlling the DFS algorithm */
			dfs_opts opts;
			
			////////////////////////////////////////////////////////////
			// Thread-local algorithm data
			////////////////////////////////////////////////////////////
			
			/** Global map of seen cobases, up to symmetry */
			cobasis_map basisOrbits;
			/** Cache of recently seen cobases */
			lru::cache<index_set, index_set_hash> cobasisCache;
			/** Lookup cobases by gram vector */
			cobasis_gram_map cobasisGramMap;
			/** representatives of each orbit (of rays) */
			coordinates_map rayOrbits;
			/** Index in the global ray orbit list this explorer is updated
			 *  to */
			uind rayUpdate;
			/** representatives of each orbit (of vertices) */
			coordinates_map vertexOrbits;
			/** Lookup vertices by gram vector */
			vertex_gram_map vertexGramMap;
		};
	
	public:
		
		////////////////////////////////////////////////////////////////
		// Public interface
		////////////////////////////////////////////////////////////////
		
		/** Set up a DFS on the given matrix, with the given permuation 
		 *  group.
		 *  @param m		The matrix to DFS on
		 *  @param lin		The set of indices that are linearities
		 *  @param g		The permutation group of the matrix
		 *  @param gram		The gram matrix for the constraints (may be 
		 *  				empty if o.gramVec == false)
		 *  @param opts		The options for this DFS (default values if 
		 *  				not provided)
		 */
		dfs(matrix& m, index_set& lin, permutation_group& g, 
			gram_matrix& gram, dfs_opts o = dfs_opts());
		
		/** Perform the DFS. Should be called single-threaded, handles 
		 *  parallelism internally.
		 *  @return did the DFS run to completion, or halt because of 
		 *  	too many bases?
		 */
		bool doDfs();
		
		////////////////////////////////////////////////////////////////
		// Query methods for after completion of doDfs()
		////////////////////////////////////////////////////////////////
		
		/** @return representatives of each of the orbits of the 
		 *  cobases */
		cobasis_map getBasisOrbits() const;
		
		/** @return the dimension of the polytope */
		ind getDimension() const;
		
		/** @return the initial cobasis for the DFS */
		index_set getInitialCobasis() const;
		
		/** @return did the DFS complete (true) or terminate due to too 
		 *  many bases (false) */
		bool isFinished() const;
		
		/** @return representatives of each of the orbits of the 
		 *  extreme rays */
		coordinates_map getRayOrbits() const;
		
		/** @return the running time of the algorithm, in 
		 *  milliseconds */
		std::clock_t getRunningTime() const;
		
		/** @return the symmetry group used in the DFS */
		permutation_group const& getSymmetryGroup() const;
		
		/** @return representatives of each of the orbits of the 
		 *  vertices */
		coordinates_map getVertexOrbits() const;
		
		/** @return the gram matrix. (This method is of primary use for 
		 *  debugging the gram vectors.) */
		gram_matrix const& getGramMat() const;
		
	private:
		
		////////////////////////////////////////////////////////////////
		// Private interface
		////////////////////////////////////////////////////////////////
		
		/** Adds a cobasis to the global map of cobasis 
		 *  representatives. Note that the vertex data pointer given 
		 *  should already be present in the global map of vertex 
		 *  representatives (added by another cobasis), or addVertex() 
		 *  should be used (which calls this method). This should not 
		 *  be called from within a critical section, as it handles its 
		 *  own synchronization.
		 *  @param cob		The cobasis to map to the vertex
		 *  @param dat		The data for this vertex
		 */
		void addCobasis(index_set const& cob, vertex_data_ptr dat);
		
		/** Adds a vertex to the global map of vertex representatives.
		 *  This will add all cobases defined for this vertex to the 
		 *  global map of cobasis representatives as well. This should 
		 *  not be called from within a critical section, as it handles 
		 *  its own synchronization.
		 *  @param dat		The data for this vertex
		 */
		void addVertex(vertex_data_ptr dat);
		
		/** Gets the gram vector for an incidence set.
		 *  @param inc		The incidence set to take the gram vector 
		 *  				for
		 *  @return the gram matrix, restricted to this incidence set 
		 *  		in row and column indices, then sorted
		 */
		gram_matrix fastGramVec(index_set inc);
		
		/** Finds the rays in the current dictionary of the explorer.
		 *  @param ex		The explorer to use to find the rays
		 */
		void getRays(explorer ex);

		/** Initializes algorithm globals */
		void initGlobals();
		
		/** Combines the cobasis and coordinate data into a vertex_data 
		 *  object
		 *  @param cob		The cobasis data
		 *  @param coords	The ray coordinates
		 *  @return the vertex data from the parameters
		 */
		vertex_data_ptr rayData(cobasis_ptr cob, vector_mpz_ptr coords);
		
		/** Combines the cobasis and coordinate data into a vertex_data 
		 *  object
		 *  @param cob		The cobasis data
		 *  @param coords	The vertex coordinates (un-rationalized)
		 *  @return the vertex data from the parameters
		 */
		vertex_data_ptr vertexData(cobasis_ptr cob, 
				vector_mpz_ptr coords);
		
		////////////////////////////////////////////////////////////////
		// Initialization-time globals
		////////////////////////////////////////////////////////////////
		
		/** Matrix to store */
		matrix& m;
		/** Index set of linearities */
		index_set& lin;
		/** Permutation group used for this DFS */
		permutation_group& g;
		/** Options for controlling the DFS algorithm */
		dfs_opts opts;
		/** Dimension of the problem */
		ind dim;
		/** number of rows in the problem */
		ind rows;
		/** matrix of pre-computed inner product representatives, for 
		 *  gram vectors */
		gram_matrix gramMat;
		
		////////////////////////////////////////////////////////////////
		// Algorithm Data
		////////////////////////////////////////////////////////////////
		
		/** Global map of seen cobases, up to symmetry */
		cobasis_list globalBasisOrbits;
		/** Temporary to store time diffs into. Will hold total running 
		 *  time on algorithm completion. */
		std::clock_t diff_time;
		/** If the basis count hit the maximum count */
		bool hitMaxBasis;
		/** The first cobasis found */
		index_set initialCobasis;
		/** representatives of each orbit (of rays) */
		coordinates_list globalRayOrbits;
		/** The true dimension of the polytope */
		ind realDim;
		/** the time at which the algorithm was started */
		std::clock_t start_time;
		/** representatives of each orbit (of vertices) */
		coordinates_list globalVertexOrbits;
				
		////////////////////////////////////////////////////////////////
		// Time-related functions and values
		////////////////////////////////////////////////////////////////
		
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