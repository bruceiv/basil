#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>

#include <gmp.h>
#include <gmpxx.h>

#include <permlib/permlib_api.h>

#include "dfs.hpp"
#include "fmt.hpp"

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
		
		/** Checks if a string begins with a given prefix
		 *  @param s		the string to search
		 *  @param pre		the prefix to match
		 *  @return true if the prefix matches, false otherwise
		 */
		bool prefixMatch(std::string const& s, char const* pre) {
			for (unsigned int i = 0; pre[i]; ++i) {
				if ( 
					i >= s.length() 	/* s is shorter than prefix */
					|| s[i] != pre[i]	/* mismatch with prefix */
				) return false;
			}
			return true;				/* match with prefix in all locations */
		}
		
		/** Allocates new matrix on heap.
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
		void genMatrixFromStream(std::istream& in) {
			
			string s = "";
			/* temporary linearity vector */
			std::vector<ind> linV(0);
			
			/* parse options up to begin line */
			while ( ! prefixMatch(s, "begin") ) {
				
				if ( prefixMatch(s, "V-representation") ) {
					/* Set vertex representation flag */
					dfsOpts_.inVRepresentation();
					if (verbose) out() << "**V-representation**" << std::endl;
				}
				
				if ( prefixMatch(s, "A-representation") ) {
					/* Set vertex representation flag */
					dfsOpts_.inARepresentation();
					if (verbose) out() << "**A-representation**" << std::endl;
				}
				
				if ( prefixMatch(s, "linearity") ) {
					/* parse linearities */
					std::istringstream read(s);
					ind k, t;
					read >> k; /* linearity count */
					linV.resize(k);
					/* read linearities */
					for (ind i = 0; i < k; ++i) { read >> linV[i]; }
				}
				
				/* get next line */
				getContentLine(in, s);
			}
			
			/* get dimension line */
			getContentLine(in, s);
			std::istringstream read(s);
			
			/* read dimensions */
			ind n, d;
			read >> n;
			read >> d;
			
			/* create new matrix and load data */
			m = boost::make_shared<matrix>(n, d);
			for (ind i = 0; i < n; i++) {
				for (ind j = 0; j < d; j++) {
					in >> m->elem(i,j);
					m->elem(i,j).canonicalize();
				}
			}
			
			/* ignore up to end line */
			getContentLine(in, s);
			while ( s != string("end") ) getContentLine(in, s);
			
			/* read linearities into index set */
			l = dfs::index_set(n+1);
			for (std::vector<ind>::iterator iter = linV.begin(); 
					iter != linV.end(); ++iter) l.set(*iter);
			
			/* NOTE if not split input, could consume the rest of it here to 
			 * parse options. */
			
		}
	
		/** Allocates new permutation group on heap and returns it.
		 *  creates permutation groups on n elements, where n is m.size1()
		 *  expects input on the stream in to be a newline-delimited list of 
		 *  permuataions, where a permutation is a comma-delimited lists of 
		 *  cycles, and a cycle is a whitespace-delimited list of elements from 
		 *  the range [1..n].
		 */
		void genPermutationGroupFromStream(std::istream& in) {
			
			ind n = m->size();
			
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
			
			g = permlib::construct(n, generators.begin(), generators.end());
		}
		
	public:
		/** Argument constructor; parses program arguments. Read code for 
		 *  option descriptions
		 *  
		 *  @param argc		The number of arguments
		 *  @param argv		The list of arguments
		 */
		opts(int argc, char** argv) 
				: dfsOpts_(), matFile(), grpFile(), outFile(), 
				splitInput(false), m(), g(), l(), verbose(true) {
			
			using namespace boost::program_options;
			
			string matFileName = "", grpFileName = "", outFileName = "";
			long printInterval = 0;
			
			options_description o("Basil options");
			o.add_options()
				("arrangement-pivot", 
					bool_switch(&dfsOpts_.aRepresentation),
					"Makes Basil pivot as if the input was an arrangement.")
				("assume-no-symmetry", 
					bool_switch(&dfsOpts_.assumesNoSymmetry),
					"Forces Basil to assume there is no symmetry in the input.")
				("show-all-dicts", 
					bool_switch(&dfsOpts_.showsAllDicts), 
					"Show all intermediate dictionaries in the search tree.")
				("no-gram-vec", 
					bool_switch(&dfsOpts_.gramVec)
						->default_value(true)->implicit_value(false),
					"Deactivate gram vector hashing (not reccommended, gram "
					"vectors are a cheap optimization)")
				("stab-search", 
					bool_switch(&dfsOpts_.stabSearch),
					"Activate cobasis stabilizer search (not reccommended, "
					"stabilizer computation costs more than it saves)")
				("print-basis",
					value<long>(&dfsOpts_.printBasis),
					"Print the number of cobases found and running time every "
					"n cobases.")
				("print-interval",
					value<long>(&printInterval),
					"Convenience for print-{basis,ray,vertex} with the given "
					"parameter. If any of the others are given, they take "
					"precedence")
				("print-new",
					bool_switch(&dfsOpts_.printNew),
					"Print the added {cobasis,vertex,ray} when printing a "
					"status message")
				("print-ray",
					value<long>(&dfsOpts_.printRay),
					"Print the number of cobases found and running time every "
					"n cobases.")
				("print-vertex",
					value<long>(&dfsOpts_.printVertex),
					"Print the number of cobases found and running time every "
					"n cobases.")
				("verbose,v",
					bool_switch(&verbose),
					"Shorthand for --print-interval=1, --print-new. Those "
					"settings, if supplied, will take precedence.")
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
			
			/* Handle verbose flag */
			if ( verbose ) {
				if ( ! printInterval ) printInterval = 1;
				dfsOpts_.doPrintNew();
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
			
			genMatrixFromStream(matIn());
			genPermutationGroupFromStream(grpIn());
		}
		
		/** Destructor. */
		~opts() {
			if ( outFile.is_open() ) outFile.close();
			if ( grpFile.is_open() ) grpFile.close();
			if ( matFile.is_open() ) matFile.close();
		}
		
		/** get the output stream */
		std::ostream& out()
			{ return outFile.is_open() ? outFile : std::cout; }
		
		/** get the DFS options */
		dfs_opts& dfsOpts() { return dfsOpts_; }
		
		/** get the problem matrix. Will read it from the matrix input stream 
		 *  the first time, caching it for later use */
		matrix& mat() {
			return *m;
		}
		
		/** get the problem permutation group. Will read it from the group 
		 *  input stream the first time, caching it for later use */
		permutation_group& grp() {
			return *g;
		}
		
		dfs::index_set& lin() {
			return l;
		}
		
		/** true if input is split over two files, false otherwise */
		bool isVerbose() { return verbose; }
		
	private:
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
		
		/** options to provide to the DFS */
		dfs_opts dfsOpts_;
		/** Matrix input stream */
		std::ifstream matFile;
		/** Group input stream */
		std::ifstream grpFile;
		/** Output stream */
		std::ofstream outFile;
		/** true if the input is split across two files [false]. Can be made 
		 *  true by specifying group file, or explicitly specifying group */
		bool splitInput;
		
		/** Pointer to the matrix for the problem */
		matrix_ptr m;
		/** Pointer to the permutation group for the problem */
		permutation_group_ptr g;
		/** linearity indices */
		dfs::index_set l;
		
		/** verbose output printing [false]. */
		bool verbose;
	}; /* class opts */
} /* namespace basil */


/** Test driver for basil.
 *  Accepts a matrix and permutation_group on standard input, generating them 
 *  using genMatrixFromStream() and genPermutationGroupFromStream().
 */
int main(int argc, char **argv) {
	using namespace basil;
	
	/* parse command line arguments */
	opts o(argc, argv);
	
	std::ostream& out = o.out();
	std::ostream& (*endl)(std::ostream&) = std::endl;
	
	if ( o.isVerbose() ) {
		//read in & print matrix
		out << "Matrix:\t" <<  fmt( o.mat(), 0 ) << endl;
		
		//read in & print permutation group
		out << "Group:\t" << fmt( o.grp(), 0 ) << endl;
	}
	
	//initialize DFS algorithm
	dfs d(o.mat(), o.lin(), o.grp(), o.dfsOpts() );
	
	if ( o.isVerbose() ) {
		//print inner product matrix (NOTE for debugging)
		out << "Inner Product Matrix:\t" << fmt( d.getInnerProdMat(), 0 ) 
				<< endl;
	}
	
	//run DFS algorithm
	if ( d.doDfs() ) {
		
		out 	<< "\nresults: " 
				<< "\n{"
				<< "\n\tdimension: " << d.getDimension()
				<< "\n\tinitial cobasis: " << fmt( d.getInitialCobasis() )
				<< "\n\tsymmetry generators: " << fmt( d.getSymmetryGroup(), 1 )
				<< "\n\tbasis orbits: " << fmt( d.getBasisOrbits(), 1 )
				<< "\n\tvertex orbits: " << fmt( d.getVertexOrbits(), 1 )
				<< "\n\tray orbits: " << fmt( d.getRayOrbits(), 1 )
				<< "\n}"
				<< "\ntotal running time: " << d.getRunningTime() << " ms"
				<< endl;
		
	} else {
		out << "ERROR: DFS terminated due to too many bases" << endl;
	}
		
	return 0;
}
