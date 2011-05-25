#include <iostream>

#include "basilCommon.hpp"
#include "dfs.hpp"
#include "groupGen.hpp"
#include "matrixGen.hpp"

using namespace std;
using namespace basil;

int main(int argc, char **argv) {
    
	//read in & print matrix
	shared_ptr<matrix> m = genMatrixfromStream(cin);
	cout << *m << endl;
	
	//read in & print permutation group
	shared_ptr<permutation_group> g = genPermutationGroupFromStream(cin, *m);
	cout << *g << endl;
	
	//initialize DFS algorithm
	dfs d;
	
	return 0;
}
