#include <iostream>

#include "basilCommon.hpp"
#include "dfs.hpp"
#include "groupGen.hpp"
#include "matrixGen.hpp"

using namespace std;
using namespace basil;

/** Simple test driver for basil.
 *  Accepts a matrix and permutation_group on standard input, generating them 
 *  using genMatrixFromStream() and genPermutationGroupFromStream().
 */
int main(int argc, char **argv) {
    
	//read in & print matrix
	shared_ptr<matrix> m(genMatrixFromStream(cin));
	
	cout << *m << endl;
	
	//read in & print permutation group
	shared_ptr<permutation_group> g(genPermutationGroupFromStream(cin, *m));
	cout << *g << endl;
	
	//initialize DFS algorithm
	dfs d(*m, *g, dfs_opts().showsAllDicts() );
	
	//run DFS algorithm
	d.doDfs();
		
	return 0;
}
