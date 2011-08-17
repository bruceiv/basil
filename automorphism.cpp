#include <algorithm>
#include <list>
#include <set>

#include <boost/make_shared.hpp>

#include <permlib/common.h>
#include <permlib/sorter/base_sorter.h>
//above included because the PermLib author didn't ...
#include <permlib/symmetric_group.h>
#include <permlib/search/partition/matrix_automorphism_search.h>

#include "automorphism.hpp"
#include "basil.hpp"
#include "gram.hpp"
#include "permUtils.hpp"

namespace basil {
	
	permutation_group_ptr compute_restricted_automorphisms(
			const gram_matrix& g) {
		
		/* typedefs */
		typedef
			permlib::SymmetricGroup<permutation>
			symmetric_group;
		
		typedef
			permlib::partition::MatrixAutomorphismSearch<
				symmetric_group, 
				permutation_transversal
				>
			matrix_automorphism_search;
		
		/* degree of the permutation groups (size of the input matrix) */
		uind n = g.dim();
		
		/* overall group to search in - the symmetric group over n */
		symmetric_group s_n(n);
		
		/* canonicalize g for permlib */
		gram_matrix c = g.permlibCanon();
		
		/* the search to perform */
		matrix_automorphism_search mas(s_n, 0);
		mas.construct(c);
		
		/* perform search */
		permutation_group_ptr ret = boost::make_shared<permutation_group>(n);
		mas.search(*ret);
		
		/* return results */
		return ret;
	}
	
	permutation_group_ptr compute_arrangement_automorphisms(
			const gram_matrix& g) {
		
		/* compute automorphisms of sign-doubled gram matrix */
		permutation_group_ptr g_d = 
				compute_restricted_automorphisms(g.doubled());
		
// std::cout << "\nSign-doubled group:";
// permutation_list sgs = strong_gen_set(*g_d);
// for (permutation_list::iterator si = sgs.begin(); si != sgs.end(); ++si)
// std::cout << " " << **si;
		/* 
		 * reverse sign-doubling process on the output 
		 */
		
		/* generators of the final permutation group */
		std::vector<permutation_ptr> gens;
		
		/* for each generator of the sign-doubled permutation group */
		permutation_list gen_list = strong_gen_set(*g_d);
		for (permutation_list::iterator iter = gen_list.begin(); 
			 iter != gen_list.end(); ++iter) {
			
			/* sign-doubled permutation */
			permutation_cycle_list p_d = cycle_list(**iter);
			
			/* list of expected cycles that are negations of seen permutations 
			 * (and thus redundant) */
			std::list<permutation_cycle> red;
			
			/* list of cycles that have been un-sign doubled */
			permutation_cycle_list p;
			
			for (permutation_cycle_list::iterator cIter = p_d.begin(); 
				 cIter != p_d.end(); ++cIter) {
				
				/* sign-doubled cycle */
				permutation_cycle c_d = *cIter;
				
				/* ignore the cycle if you can find it in the expected list */
				bool ignore = false;
				for (std::list<permutation_cycle>::iterator lIter = red.begin();
					 lIter != red.end(); ++lIter) {
					if ( std::equal(c_d.begin(), c_d.end(), lIter->begin()) ) {
						ignore = true; break;
					}
				}
				if (ignore) continue;
				
				/* negation of cycle */
				permutation_cycle c_n(c_d.size());
				for (uind i = 0; i < c_d.size(); ++i) {
					/* toggling the low order bit switches a row for its 
					 * negation in this sign-doubling scheme */
					c_n[i] = c_d[i] ^ 0x1;
				}
				
				permutation_cycle::iterator nIter = 
						std::find(c_d.begin(), c_d.end(), c_n.front());
				if ( nIter == c_d.end() ) {
					/* if the negated cycle is disjoint, add it to the 
					 * redundant list */
					red.push_back(c_n);
				} else if ( c_d.size() == 2 ) {
					/* if the cycle is of the form (x,-x), ignore it */
					continue;
				}
				
				/* cycle in the original space */
				permutation_cycle c;
				for (permutation_cycle::iterator eIter = c_d.begin(); 
					 eIter != nIter; ++eIter) {
					/* removing the low order bit takes a row to it's original 
					 * representation in this sign-doubling scheme */
					c.push_back( (*eIter) >> 1 );
				}
				
				/* add cycle to the permutation */
				p.push_back(c);
			}
			
			/* add permutation to generator list */
			gens.push_back(boost::make_shared<permutation>(perm(g.dim(), p)));
		}
		
// std::cout << "\nReconstituted group:";
// for (std::vector<permutation_ptr>::iterator gi = gens.begin(); gi != gens.end(); ++gi)
// std::cout << " " << **gi;
		return permlib::construct(g.dim(), gens.begin(), gens.end());
	}
	
} /* namespace basil */
