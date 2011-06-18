#include <iostream>
#include <fstream>

#include "basilCommon.hpp"
#include "dfs.hpp"
#include "groupGen.hpp"
#include "matrixGen.hpp"

using namespace basil;

/** Prints a representation of its cobasis (as a set of indices).
 *  NOTE this redefines the standard output operator for an index set, and all 
 *  code included afterward will use this operator. Perhaps this should be 
 *  refactored ...
 */
std::ostream& operator<< (std::ostream& o, lrs::index_set const& s) {
	bool isFirst = true;
	o << "{";
	for (lrs::index_set_iter it = lrs::begin(s); 
			it != lrs::end(s); 
			++it) {
		if (isFirst) isFirst = false; else o << ", ";
		o << *it;
	}
	o << "}";
	return o;
}

/** Prints a list of cobases. */
std::ostream& operator<< (std::ostream& o, dfs::cobasis_invariants_list l) {
	bool isFirst = true;
	o << "{";
	for (dfs::cobasis_invariants_list::const_iterator it = l.begin();
			it != l.end(); ++it) {
		if (isFirst) isFirst = false; else o << ", ";
		o << (**it).cob;
	}
	o << "}";
	return o;
}

/** prints a representation of a list of permutations */
std::ostream& operator<< (std::ostream& o, permutation_list const& l) {
	bool isFirst = true;
	o << "{";
	for (permutation_list::const_iterator it = l.begin();
			it != l.end(); ++it) {
		if (isFirst) isFirst = false; else o << ", ";
		o << **it;
	}
	o << "}";
	return o;
}

/** prints a representation of a list of vertices */
std::ostream& operator<< (std::ostream& o, dfs::vertex_rep_list const& l) {
	bool isFirst = true;
	o << "{";
	for (dfs::vertex_rep_list::const_iterator it = l.begin();
			it != l.end(); ++it) {
		if (isFirst) isFirst = false; else o << ", ";
		o << (**it).coords;
	}
	o << "}";
	return o;
}

namespace basil {
	/** Contains runtime custimizations for basil
	 */
	class opts {
	public:
		/** Default constructor. Equivalent to calling iterator constructor with an 
		*  empty iterator */
		opts() : dfsOpts_() {
			/* set default options */
			dfsOpts_.showAllDicts();
		}
		
		/** Iterator constructor; parses program arguments.
		*  Arguments are expected in the following format:
		*  	./basil [ matIn [grpIn] ]
		*  
		*  matIn: the name of the matrix input file [default: standard input]
		*  grpIn: the name of the permutation group input file [default: matIn]
		*  
		*  @param begin	The beginning iterator for the list of arguments
		*  @param end		The end iterator for the list of arguments
		*  @param Iter		the type of the options iterator
		*/
		template <typename Iter>
		opts(Iter begin, Iter end) : dfsOpts_() {
			/* consume program name */
			if (begin != end) ++begin;
			
			/* set default options */
			dfsOpts_.showAllDicts();
			
			if (begin != end && str_equals(*begin, "--assume-no-symmetry")) {
				dfsOpts_.assumeNoSymmetry();
				++begin;
			}
			
			/* parse matrix file name */
			if (begin != end) {
				matFile.open(*begin);
				++begin;
			}
			
			/* parse group file name */
			if (begin != end) {
				grpFile.open(*begin);
				++begin;
			}
		}
		
		/** Destructor. */
		~opts() {
			if ( grpFile.is_open() ) grpFile.close();
			if ( matFile.is_open() ) matFile.close();
		}
		
		/** get the input stream for the matrix */
		std::istream& matIn() 
			{ return matFile.is_open() ? matFile : std::cin; }
		
		std::istream& grpIn()
			{ return grpFile.is_open() ? grpFile : matIn(); }
		
		dfs_opts& dfsOpts() { return dfsOpts_; }
		
	private:
		
		/** Checks two objects for equality by casting them to std::string. */
		template<typename T1, typename T2>
		bool str_equals(T1 a, T2 b) {
			return string(a) == string(b);
		}
		
		/** options to provide to the DFS */
		dfs_opts dfsOpts_;
		/** Matrix input stream */
		std::ifstream matFile;
		/** Group input stream */
		std::ifstream grpFile;
	}; /* class opts */
} /* namespace basil */


/** Simple test driver for basil.
 *  Accepts a matrix and permutation_group on standard input, generating them 
 *  using genMatrixFromStream() and genPermutationGroupFromStream().
 */
int main(int argc, char **argv) {
	/* parse command line arguments */
	opts o(argv, argv+argc);
	
	std::ostream& cout = std::cout;
	std::ostream& (*endl)(std::ostream&) = std::endl;
	
	//read in & print matrix
	matrix_ptr m(genMatrixFromStream( o.matIn() ));
	cout << *m << endl;
	
	//read in & print permutation group
	permutation_group_ptr g(genPermutationGroupFromStream(o.grpIn(), *m));
	cout << *g << endl;
	
	//initialize DFS algorithm NOTE debug mode
	dfs d(*m, *g, o.dfsOpts() );
	
	//run DFS algorithm
	if ( d.doDfs() ) {
		
		cout 	<< "\nBEGIN RESULTS" 
				<< "\n{"
				<< "\n\tdimension:" << d.getDimension()
				<< "\n\tinitial cobasis: " << d.getInitialCobasis() 
				<< "\n\tsymmetry generators: " << d.getSymmetryGroup().S 
				<< "\n\tbasis orbits: " << d.getBasisOrbits()
				<< "\n\tvertex orbits: " << d.getVertexOrbits()
				<< "\n\tray orbits: " << d.getRayOrbits()
				<< "\n}"
				<< "\nEND RESULTS" 
				<< endl;
		
	} else {
		cout << "ERROR: DFS terminated due to too many bases" << endl;
	}
		
	return 0;
}
