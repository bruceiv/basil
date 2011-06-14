#include <iostream>

#include "basilCommon.hpp"
#include "dfs.hpp"
#include "groupGen.hpp"
#include "matrixGen.hpp"

using namespace basil;

/** Prints a representation of its cobasis (as a set of indices). */
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


/** Simple test driver for basil.
 *  Accepts a matrix and permutation_group on standard input, generating them 
 *  using genMatrixFromStream() and genPermutationGroupFromStream().
 */
int main(int argc, char **argv) {
	std::istream& cin = std::cin;
	std::ostream& cout = std::cout;
	std::ostream& (*endl)(std::ostream&) = std::endl;
	
	//read in & print matrix
	matrix_ptr m(genMatrixFromStream(cin));
	cout << *m << endl;
	
	//read in & print permutation group
	permutation_group_ptr g(genPermutationGroupFromStream(cin, *m));
	cout << *g << endl;
	
	//initialize DFS algorithm NOTE debug mode, no PermLib
	dfs d(*m, *g, dfs_opts().showAllDicts() );
	
	//run DFS algorithm
	if ( d.doDfs() ) {
				
		cout << "\nBEGIN RESULTS" 
				<< "\n{"
				<< "\n\tdimension:" << d.getDimension()
				<< "\n\tinitial cobasis: " << d.getInitialCobasis() 
				<< "\n\tsymmetry generators: " << d.getSymmetryGroup().S 
				<< "\n\tbasis orbits #: " << d.getBasisOrbits().size()
				<< "\n\tvertex orbits #: " << d.getVertexOrbits().size()
				<< "\n\tray orbits #: " << d.getRayOrbits().size()
				<< "\n}"
				<< "\nEND RESULTS" << endl;
	} else {
		cout << "ERROR: DFS terminated due to too many bases";
	}
		
	return 0;
}
