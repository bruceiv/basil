#ifndef _GROUP_GEN_HPP_
#define _GROUP_GEN_HPP_

#include <istream>
#include <vector>

#include <permlib/permlib_api.h>

#include "basilCommon.hpp"

namespace basil {
	
	/** Allocates new permutation group on heap and returns it.
	 *  creates permutation groups on n elements, where n is m.size1()
	 *  expects input on the stream in to be a newline-delimited list of 
	 * permuataions, where a permutation is a comma-delimited lists of cycles, 
	 * and a cycle is a whitespace-delimited list of elements from the range 
	 * [1..n]
	 */
	shared_ptr<permutation_group> genPermutationGroupFromStream(
			std::istream& in, const matrix& m) {
		
		using namespace std;
		
		ind n = m.size1();
		
		//vector< shared_ptr<permutation> > generators;
		vector< shared_ptr<permutation> > generators;
		
		//read in generators
		string s;
		shared_ptr<permutation> p;
		while (getline(in, s)) {
			p.reset(new permutation(n, s));
			generators.push_back(p);
		}
		
		return permlib::construct(n, generators.begin(), generators.end());
		
	}
	
} /* namespace basil */

#endif /* _GROUP_GEN_HPP_ */