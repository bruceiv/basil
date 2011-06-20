#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

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
		opts() : dfsOpts_() {}
		
		/** Iterator constructor; parses program arguments.
		*  Arguments are expected in the following format:
		*  	./basil [ matIn [grpIn] ]
		*  
		*  matIn: the name of the matrix input file [default: standard input]
		*  grpIn: the name of the permutation group input file [default: matIn]
		*  
		*  @param argc		The number of arguments
		*  @param argv		The list of arguments
		*/
		opts(int argc, char** argv) : dfsOpts_() {
			using namespace boost::program_options;
			
			string matFileName = "", grpFileName = "", outFileName = "";
			
			options_description o("Basil options");
			o.add_options()
				("assume-no-symmetry", 
					bool_switch(&dfsOpts_.assumesNoSymmetry),
					"Forces Basil to assume there is no symmetry in the input.")
				("show-all-dicts", 
					bool_switch(&dfsOpts_.showsAllDicts), 
					"Show all intermediate dictionaries in the search tree.")
				("input-file,i",
					value<string>(&matFileName),
					"Input file name. "
					"Alias for matrix-file.")
				("matrix-file,m",
					value<string>(&matFileName),
					"Matrix file name. "
					"Standard input if none supplied; may also be supplied as "
					"first positional argument.")
				("group-file,g",
					value<string>(&grpFileName),
					"Group file name. "
					"Matrix input stream if none supplied; may also be "
					"supplied as second positional argument.")
				("output-file,o",
					value<string>(&outFileName),
					"Output file name. "
					"Standard output if none supplied; may also be supplied as "
					"third positional argument.")
				;
			positional_options_description p;
 			p.add("input-file", 1).add("group-file", 2).add("output-file", 3);
			
			variables_map v;
			store(command_line_parser(argc, argv)
					.options(o).positional(p).allow_unregistered().run(), v);
			notify(v);
			
			/* Open I/O files, if supplied */
			if ( ! matFileName.empty() ) matFile.open(matFileName.c_str());
			if ( ! grpFileName.empty() ) grpFile.open(grpFileName.c_str());
			if ( ! outFileName.empty() ) outFile.open(outFileName.c_str());
		}
		
		/** Destructor. */
		~opts() {
			if ( outFile.is_open() ) outFile.close();
			if ( grpFile.is_open() ) grpFile.close();
			if ( matFile.is_open() ) matFile.close();
		}
		
		/** get the input stream for the matrix */
		std::istream& matIn() 
			{ return matFile.is_open() ? matFile : std::cin; }
		
		std::istream& grpIn()
			{ return grpFile.is_open() ? grpFile : matIn(); }
		
		std::ostream& out()
			{ return outFile.is_open() ? outFile : std::cout; }
		
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
		/** Output stream */
		std::ofstream outFile;
	}; /* class opts */
} /* namespace basil */


/** Simple test driver for basil.
 *  Accepts a matrix and permutation_group on standard input, generating them 
 *  using genMatrixFromStream() and genPermutationGroupFromStream().
 */
int main(int argc, char **argv) {
	/* parse command line arguments */
	opts o(argc, argv);
	
	std::ostream& out = o.out();
	std::ostream& (*endl)(std::ostream&) = std::endl;
	
	//read in & print matrix
	matrix_ptr m(genMatrixFromStream( o.matIn() ));
	out << *m << endl;
	
	//read in & print permutation group
	permutation_group_ptr g(genPermutationGroupFromStream(o.grpIn(), *m));
	out << *g << endl;
	
	//initialize DFS algorithm NOTE debug mode
	dfs d(*m, *g, o.dfsOpts() );
	
	//run DFS algorithm
	if ( d.doDfs() ) {
		
		out 	<< "\nBEGIN RESULTS" 
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
		out << "ERROR: DFS terminated due to too many bases" << endl;
	}
		
	return 0;
}
