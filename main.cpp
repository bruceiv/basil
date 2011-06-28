#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/program_options.hpp>

#include <gmpxx.h>

#include <permlib/permlib_api.h>

#include "dfs.hpp"


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
std::ostream& operator<< (std::ostream& o, dfs::cobasis_map const& m) {
	/* sort the set of cobases */
	std::set<dfs::index_set> s;
	for (dfs::cobasis_map::const_iterator it = m.begin(); it != m.end(); ++it) {
		s.insert(it->first);
	}
	
	/* print the sorted set */
	bool isFirst = true;
	o << "{";
	for (std::set<dfs::index_set>::const_iterator it = s.begin(); 
			it != s.end(); ++it) {
		if (isFirst) isFirst = false; else o << ", ";
		o << *it;
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
std::ostream& operator<< (std::ostream& o, dfs::coordinates_map const& m) {
	/* sort the set of coordinates */
	std::set<dfs::coordinates> s;
	for (dfs::coordinates_map::const_iterator it = m.begin(); it != m.end(); 
			++it) {
		s.insert(it->first);
	}
	
	/* print the sorted set */
	bool isFirst = true;
	o << "{";
	for (std::set<dfs::coordinates>::const_iterator it = s.begin(); 
			it != s.end(); ++it) {
		if (isFirst) isFirst = false; else o << ", ";
		o << *it;
	}
	o << "}";
	return o;
}

namespace basil {
	/** Contains runtime custimizations for basil
	 */
	class opts {
	private:
		/* utility methods */
		
		/** Gets a line that isn't a comment (begins with '*' or '#').
		 *  @param in		the input stream to read the line from
		 *  @param s		the string to store the line in
		 *  @return the input stream
		 */
		std::istream& getContentLine(std::istream& in, string& s) {
			do {
				std::getline(in, s);
			} while (s.empty() || s[0] == '*' || s[0] == '#');
			
			return in;
		}
		
		/** Allocates new matrix on heap and returns it.
		 *  Expects input in the following format (to match lrs):
		 *  
		 *  [name]
		 *  [V-representation]
		 *  [\<lrs options\>]
		 *  begin
		 *  \<n\> \<d\> rational
		 *  \< n * d whitespace-delimited data values \>
		 *  end
		 *  
		 *  where [] denotes an optional value, \<\> denotes a variable, and any 
		 *  line beginning with '#' or '*' is ignored as a comment line.
		 */
		matrix_ptr genMatrixFromStream(std::istream& in) {
			
			string s = "";
			
			/* parse options up to begin line */
			while ( s != string("begin") ) {
				
				if ( s == string("V-representation") ) {
					/* Set vertex representation flag */
					dfsOpts_.inVRepresentation();
					out() << "**V-representation**" << std::endl;
				}
				
				/* get next line */
				getContentLine(in, s);
			}
			
			/* get dimension line */
			getContentLine(in, s);
			std::istringstream lin(s);
			
			/* read dimensions */
			ind n, d;
			lin >> n;
			lin >> d;
			
			/* create new matrix and load data */
			matrix_ptr m(new matrix(n, d));
			mpq_class t;
			for (ind i = 0; i < n; i++) {
				for (ind j = 0; j < d; j++) {
					in >> t;
					(*m)[i][j] = t;
				}
			}
			
			/* ignore up to end line */
			getContentLine(in, s);
			while ( s != string("end") ) getContentLine(in, s);
			
			/* NOTE if not split input, could consume the rest of it here to 
			 * parse options. */
			
			return m;
		}
	
		/** Allocates new permutation group on heap and returns it.
		 *  creates permutation groups on n elements, where n is m.size1()
		 *  expects input on the stream in to be a newline-delimited list of 
		 *  permuataions, where a permutation is a comma-delimited lists of 
		 *  cycles, and a cycle is a whitespace-delimited list of elements from 
		 *  the range [1..n].
		 */
		permutation_group_ptr genPermutationGroupFromStream(std::istream& in, 
															const matrix& m) {
			
			ind n = m.n();
			
			std::vector<permutation_ptr> generators;
			
			//read in generators
			string s = "";
			permutation_ptr p;
			while ( s != string("symmetry begin") ) getContentLine(in, s);
			
			std::getline(in, s);
			while ( s != string("symmetry end") ) {
				p.reset(new permutation(n, s));
				generators.push_back(p);
				
				std::getline(in, s);
			}
			
			return permlib::construct(n, generators.begin(), generators.end());
			
		}
		
	public:
		/** Argument constructor; parses program arguments. Read code for 
		 *  option descriptions
		 *  
		 *  @param argc		The number of arguments
		 *  @param argv		The list of arguments
		 */
		opts(int argc, char** argv) 
				: dfsOpts_(), matFile(), grpFile(), outFile(), m(), g(), 
				splitInput(false) {
			
			using namespace boost::program_options;
			
			string matFileName = "", grpFileName = "", outFileName = "";
			long printInterval = 0;
			
			options_description o("Basil options");
			o.add_options()
				("assume-no-symmetry", 
					bool_switch(&dfsOpts_.assumesNoSymmetry),
					"Forces Basil to assume there is no symmetry in the input.")
				("show-all-dicts", 
					bool_switch(&dfsOpts_.showsAllDicts), 
					"Show all intermediate dictionaries in the search tree.")
				("print-basis",
					value<long>(&dfsOpts_.printBasis),
					"Print the number of cobases found and running time every "
					"n cobases.")
				("print-interval",
					value<long>(&printInterval),
					"Convenience for print-{basis,ray,vertex} with the given "
					"parameter. If any of the others are given, they take "
					"precedence")
				("print-ray",
					value<long>(&dfsOpts_.printRay),
					"Print the number of cobases found and running time every "
					"n cobases.")
				("print-vertex",
					value<long>(&dfsOpts_.printVertex),
					"Print the number of cobases found and running time every "
					"n cobases.")
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
			
			/* mark input as split if group or matrix file explicitly set */
			if (v.count("group-file") || v.count("matrix-file")) {
				splitInput = true;
			}
			
			/* Open I/O files, if supplied */
			if ( ! matFileName.empty() ) matFile.open(matFileName.c_str());
			if ( ! grpFileName.empty() ) grpFile.open(grpFileName.c_str());
			if ( ! outFileName.empty() ) {
				outFile.open(outFileName.c_str());
				dfsOpts_.withOutput(outFile);
			}
			
			/* Handle print-interval overloading */
			if ( printInterval ) {
				dfsOpts_.printAt(printInterval);
				
				if (v.count("print-basis")) 
					dfsOpts_.printBasisAt(v["print-basis"].as<long>());
				if (v.count("print-ray")) 
					dfsOpts_.printRayAt(v["print-ray"].as<long>());
				if (v.count("print-vertex")) 
					dfsOpts_.printVertexAt(v["print-vertex"].as<long>());
			}
		}
		
		/** Destructor. */
		~opts() {
			if ( outFile.is_open() ) outFile.close();
			if ( grpFile.is_open() ) grpFile.close();
			if ( matFile.is_open() ) matFile.close();
		}
		
		/** true if input is split over two files, false otherwise */
		bool isSplitInput() { return splitInput; }
		
		/** get the input stream for the matrix */
		std::istream& matIn() 
			{ return matFile.is_open() ? matFile : std::cin; }
		
		/** get the input stream for the permuation group */
		std::istream& grpIn() {
			return grpFile.is_open() ? 
					grpFile : 
					splitInput ? std::cin : matIn()
			;
		}
		
		/** get the output stream */
		std::ostream& out()
			{ return outFile.is_open() ? outFile : std::cout; }
		
		/** get the DFS options */
		dfs_opts& dfsOpts() { return dfsOpts_; }
		
		/** get the problem matrix. Will read it from the matrix input stream 
		 *  the first time, caching it for later use */
		matrix& mat() {
			if (!m) m = genMatrixFromStream(matIn());
			return *m;
		}
		
		/** get the problem permutation group. Will read it from the group 
		 *  input stream the first time, caching it for later use */
		permutation_group& grp() {
			if (!g) g = genPermutationGroupFromStream(grpIn(), mat());
			return *g;
		}
		
	private:
		/** options to provide to the DFS */
		dfs_opts dfsOpts_;
		/** Matrix input stream */
		std::ifstream matFile;
		/** Group input stream */
		std::ifstream grpFile;
		/** Output stream */
		std::ofstream outFile;
		
		/** Pointer to the matrix for the problem */
		matrix_ptr m;
		/** Pointer to the permutation group for the problem */
		permutation_group_ptr g;
		
		/** true if the input is split across two files [false]. Can be made 
		 *  true by specifying group file, or explicitly specifying group */
		bool splitInput;
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
	out << "Matrix:\t" <<  o.mat() << endl;
	
	//read in & print permutation group
	out << "Group:\t" << o.grp().S << endl;
	
	//initialize DFS algorithm NOTE debug mode
	dfs d(o.mat(), o.grp(), o.dfsOpts() );
	
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
				<< "\ntotal running time: " << d.getRunningTime() << " ms"
				<< "\nEND RESULTS" 
				<< endl;
		
	} else {
		out << "ERROR: DFS terminated due to too many bases" << endl;
	}
		
	return 0;
}
