/** Main driver for Basil.
 *
 *  @author Aaron Moss
 */

#include <cstdlib>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>

#include <gmp.h>
#include <gmpxx.h>

#ifdef BAS_MT
#include <omp.h>
#endif /* BAS_MT */

#include "automorphism.hpp"
#include "basil.hpp"

#ifdef BAS_MT
#include "dfsp.hpp"
#else
#include "dfs.hpp"
#endif /* BAS_MT */

#include "fmt.hpp"
#include "gram.hpp"
#include "metric.hpp"
#include "parse.hpp"

#include "permlib/permlib_api.h"

namespace basil {
	/** Switch between sequential and parallel DFS */
#ifdef BAS_MT
	typedef dfsp		dfs_t;
	typedef dfsp_opts	dfs_opts_t;
#else
	typedef dfs			dfs_t;
	typedef dfs_opts	dfs_opts_t;
#endif /* BAS_MT */

	/** Contains runtime custimizations for basil
	 */
	class opts {
	public:
		/** Argument constructor; parses program arguments. Read code for 
		 *  option descriptions
		 *  
		 *  @param argc		The number of arguments
		 *  @param argv		The list of arguments
		 */
		opts(int argc, char** argv) 
				: dfsOpts_(), matFile(), grpFile(), outFile(), 
				groupOverride(false), p(), verbose(false), 
				doNormalize(true), euclideanGram(true), qGram(false),
				doOrthoAugment(false), preprocessor(false), genSymmetry(false) {
			
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
				("generate-symmetry",
					bool_switch(&genSymmetry),
					"Forces Basil to generate a new symmetry group")
				("show-all-dicts", 
					bool_switch(&dfsOpts_.showsAllDicts), 
					"Show all intermediate dictionaries in the search tree.")
				("no-gram",
					bool_switch(&dfsOpts_.gramVec)
						->default_value(true)->implicit_value(false),
					"Deactivate Gram matrix hashing (Gram matrix are a cheap "
					"optimization, but their calculation may be expensive)")
				("euclidean-gram",
					bool_switch(&euclideanGram),
					"Use Euclidean metric for Gram matrix generation")
				("q-gram",
					bool_switch(&qGram),
					"Use Q-matrix metric for Gram matrix generation")
				("ortho-augment",
					bool_switch(&doOrthoAugment),
					"Use orthogonal augmentation to bring matrix to full rank "
					"for Q-matrix generation")
				("no-norm",
					bool_switch(&doNormalize)
						->default_value(true)->implicit_value(false),
					"If it is known that all symmetric vectors have the same "
					"norm, this flag will save expensive normalization "
					"calculations in gram matrix construction; use may result "
					"in duplicate orbit representatives.")
				("debug-gram",
					bool_switch(&dfsOpts_.debugGram),
					"Print gram vectors for vertices/rays/cobases that are "
					"printed")
				("stab-search", 
					bool_switch(&dfsOpts_.stabSearch),
					"Activate cobasis stabilizer search (not reccommended, "
					"stabilizer computation costs more than it saves)")
#ifdef BAS_MT
				("no-local-stack",
					bool_switch(&dfsOpts_.usesLocalStack)
						->default_value(true)->implicit_value(false),
					"Deactivate thread-local work stacks; This reduces memory "
					"usage at the cost of increased execution time.")
#endif
				("print-basis",
					value<long>(&dfsOpts_.printBasis),
					"Print the number of cobases found and running time every "
					"n cobases.")
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
				("print-each",
					value<long>(&printInterval),
					"Convenience for print-{basis,ray,vertex} with the given "
					"parameter. If any of the others are given, they take "
					"precedence")
				("print-trace",
					bool_switch(&dfsOpts_.printTrace),
					"Print the full trace of the DFS (warning: very verbose).")
				("preprocess,p", 
					bool_switch(&preprocessor), 
					"Do not DFS, simply do preprocessing work, and print "
					"normalized input file to output stream")
				("verbose,v",
					bool_switch(&verbose),
					"Shorthand for --print-interval=1, --print-new. Those "
					"settings, if supplied, will take precedence.")
				("input-file,i",
					value<string>(&matFileName),
					"Input file name. Standard input if none supplied; may "
					"also be supplied as first positional argument.")
				("matrix-file,m",
					value<string>(&matFileName),
					"Matrix file name. Alias for --input-file.")
				("group-file,g",
					value<string>(&grpFileName),
					"File to read symmetry group from - overrides any supplied "
					"in the input file.")
				("output-file,o",
					value<string>(&outFileName),
					"Output file name. Standard output if none supplied; may "
					"also be supplied as second positional argument.")
				("help,h",
					"Produce help message.")
				;
			positional_options_description p;
 			p.add("input-file", 1).add("output-file", 1);
			
			variables_map v;
			store(command_line_parser(argc, argv)
					.options(o).positional(p).allow_unregistered().run(), v);
			notify(v);
			
			if ( v.count("help") ) {
				/* print usage information and exit */
				std::cout << o << std::endl;
				exit(EXIT_FAILURE);
			}
			
			/* mark input as split if group file explicitly set */
			if ( v.count("group-file") ) {
				groupOverride = true;
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
		}
		
		/** Destructor. */
		~opts() {
			if ( outFile.is_open() ) outFile.close();
			if ( grpFile.is_open() ) grpFile.close();
			if ( matFile.is_open() ) matFile.close();
		}
		
		/** Parses the input. */
		void parse() {
			/* initial parse */
			p = basil::parse(matIn());
			
			/* set representation */
			if ( ! dfsOpts_.aRepresentation ) {
				if ( p->rep == arrangement ) dfsOpts_.inARepresentation();
				else if ( p->rep == vertex ) dfsOpts_.inVRepresentation();
			}
			
			/* get gram matrix */
			bool aRep = dfsOpts_.aRepresentation;


			if ( ! dfsOpts_.gramVec ) {
				/* do not use a Gram vector */
				p->gs = gram_omitted;
			} else if ( qGram && doOrthoAugment ) {
				/* construct a new Q-matrix-based Gram matrix, with orthogonal
				 * augmentation */
				p->gs = gram_augment;
			} else if ( euclideanGram && ! doNormalize ) {
				/* construct a new Gram matrix based on the Euclidean metric,
				 * with inner products normalized to the same length */
				p->gs = gram_inexact;
			} else if ( qGram ) {
				/* construct a new Q-matrix-based Gram matrix, without
				 * orthogonal augmentation */
				p->gs = gram_q;
			} else if ( euclideanGram ) {
				/* construct a new Gram matrix based on the Euclidean metric,
				 * without inner product normalization */
				p->gs = gram_euclid;
			} else if ( doOrthoAugment ) {
				/* construct a new Q-matrix-based Gram matrix, with orthogonal
				 * augmentation */
				p->gs = gram_augment;
			} else if ( ! doNormalize ) {
				/* construct a new Gram matrix based on the Euclidean metric,
				 * with inner products normalized to the same length */
				p->gs = gram_inexact;
			}

			matrix Qinv; matrix P; matrix_mpr N;
			switch( p->gs ) {
			case gram_omitted:
				p->gm = boost::make_shared<gram_matrix>(0);
				break;
			case gram_augment:
				Qinv = invQMat(orthoAugment(*p->m, !aRep));
				p->gm = boost::make_shared<gram_matrix>(constructGram(
						transformedInnerProdMat(*p->m, Qinv)));
				p->gs = gram_provided;
				break;
			case gram_q:
				Qinv = invQMat(*p->m);
				p->gm = boost::make_shared<gram_matrix>(constructGram(
						transformedInnerProdMat(*p->m, Qinv)));
				p->gs = gram_provided;
				break;
			case gram_euclid:
				N = normedInnerProdMat(*p->m);
				p->gm = boost::make_shared<gram_matrix>(constructGram(N));
				p->gs = gram_provided;
				break;
			case gram_inexact:
				P = innerProdMat(*p->m);
				p->gm = boost::make_shared<gram_matrix>(constructGram(P));
				p->gs = gram_provided;
				break;
			case gram_provided: /* do nothing */ break;
			}
			
			/* get symmetry group */
			if ( groupOverride ) {
				std::istream& in = grpIn();
				string s; std::getline(in, s);
				while ( s != string("symmetry begin") ) std::getline(in, s);
				p->g = parsePermutationGroup(in, p->m->size());
				p->ss = sym_provided;
			}
			if ( genSymmetry 
					|| !( p->ss == sym_provided 
						|| dfsOpts_.assumesNoSymmetry ) ) {
				if ( ! p->gs == gram_provided ) {
					matrix_mpr N = normedInnerProdMat(*p->m);
					p->gm = boost::make_shared<gram_matrix>(constructGram(N));
				}
				if ( aRep ) {
					p->g = compute_arrangement_automorphisms(*p->gm);
				} else {
					p->g = compute_restricted_automorphisms(*p->gm);
				}
				p->ss = sym_provided;
			}
			
		}
		
		/** Prints the input to the output stream, in a manner consistent with 
		 *  its input format. This is useful to perform preprocessing on input 
		 *  files.
		 */
		void print() {
			if ( ! p ) parse();
			out() << *p;
		}
		
		/** get the output stream */
		std::ostream& out()
			{ return outFile.is_open() ? outFile : std::cout; }
		
		/** get the DFS options */
		dfs_opts_t& dfsOpts() { return dfsOpts_; }
		
		/** get the problem matrix. - may call parse() if it has yet to be 
		 *  called */
		matrix& mat() {
			if ( ! p ) parse();
			return *p->m;
		}
		
		/** get the problem permutation group. - may call parse() if it has 
		 *  yet to be called */
		permutation_group& grp() {
			if ( ! p ) parse();
			return *p->g;
		}
		
		/** get the problem linearity set - may call parse() if it has yet to 
		 *  be called */
		index_set& lin() {
			if ( ! p ) parse();
			return *p->l;
		}
		
		/** gets the gram matrix - may call parse() if it has yet to be 
		 *  called */
		gram_matrix& gram() {
			if ( ! p ) parse();
			return *p->gm;
		}
		
		/** true if input is split over two files, false otherwise */
		bool isVerbose() { return verbose; }
		
		/** true if the options specify to preprocess the input only */
		bool isPreprocessor() { return preprocessor; }
		
	private:
		/** get the input stream for the matrix */
		std::istream& matIn() 
			{ return matFile.is_open() ? matFile : std::cin; }
		
		/** get the input stream for the permuation group */
		std::istream& grpIn() 
			{ return grpFile.is_open() ? grpFile : matIn(); }
		
		/** options to provide to the DFS */
		dfs_opts_t dfsOpts_;
		/** Matrix input stream */
		std::ifstream matFile;
		/** Group input stream */
		std::ifstream grpFile;
		/** Output stream */
		std::ofstream outFile;
		/** true if the group is given in its own file */
		bool groupOverride;
		
		/** Pointer to the parse results for the problem */
		parse_results_ptr p;
		
		/** verbose output printing [false]. */
		bool verbose;

		/** normalization calculation used in gram matrix construction 
		 *  [true] */
		bool doNormalize;
		/** Use Euclidean Gram matrix calculation [true] */
		bool euclideanGram;
		/** Use Q-matrix Gram matrix calculation [false] */
		bool qGram;
		/** Use orthagonal augmentation for Q-matrix Gram [false] */
		bool doOrthoAugment;

		/** only do pre-processing steps, rather than full calculation 
		 *  [false] */
		bool preprocessor;
		/** always generate new symmetry group [false] */
		bool genSymmetry;
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
	o.parse();
	
	std::ostream& out = o.out();
	std::ostream& (*endl)(std::ostream&) = std::endl;
	
	if ( o.isPreprocessor() ) {
		o.print();
		exit(EXIT_SUCCESS);
	}
	
	if ( o.isVerbose() ) {
		//print matrix
		out << "Matrix:\t" <<  fmt( o.mat(), 0 ) << endl;
		
		//print permutation group
		out << "Group:\t" << fmt( o.grp(), 0 ) << endl;
		
		//print gram matrix (if using gram matrix)
		if ( o.dfsOpts().gramVec ) out << "Gram Matrix:\t" << o.gram() << endl;
	}
	
	//initialize DFS algorithm
	dfs_t d(o.mat(), o.lin(), o.grp(), o.gram(), o.dfsOpts() );
	
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
				<< "\n}";
#ifdef BAS_MT
		#pragma omp parallel
		{
		#pragma omp master
		out 	<< "\nnumber threads: " << omp_get_num_threads();
		} /* omp parallel */
#endif /* BAS_MT */
		out 	<< "\ntotal running time: " << d.getRunningTime() << " ms";
#ifdef BAS_WALLTIME
		out		<< "\nwall time: " << d.getWallTime() << " ms";
#endif /* BAS_WALLTIME */
		out		<< endl;
		
		exit(EXIT_SUCCESS);
	} else {
		out << "ERROR: DFS terminated due to too many bases" << endl;
		exit(EXIT_FAILURE);
	}
}
