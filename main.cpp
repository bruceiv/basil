#include <iostream>

#include "basilCommon.hpp"
#include "dfs.hpp"
#include "groupGen.hpp"
#include "matrixGen.hpp"

using namespace basil;

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
	dfs d(*m, *g, dfs_opts().showAllDicts().assumeNoSymmetry() );
	
	//run DFS algorithm
	dfs::results r = d.doDfs();
	cout << "BEGIN RESULTS\n" << r << "END RESULTS" << endl;
		
	return 0;
}
