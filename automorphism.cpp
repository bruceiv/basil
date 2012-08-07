/** Implements matrix automorphism calculations from automorphism.hpp.
 *
 *  @author Aaron Moss
 */

#include <algorithm>
#include <list>
#include <set>

#include <boost/make_shared.hpp>

#include "automorphism.hpp"
#include "basil.hpp"
#include "gram.hpp"
#include "perm_utils.hpp"

#include "permlib/common.h"
#include "permlib/sorter/base_sorter.h"
//above included because the PermLib author didn't ...
#include "permlib/symmetric_group.h"
#include "permlib/search/partition/matrix_automorphism_search.h"

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
			
			/* un-sign-doubled permutation */
			permutation_cycle_list p;
			
			/* Valid sign-doubled permutations have the property that they only
			 * map a "postive" value x to a "negative" value -y where x == y */
			bool valid = true;

			for (permutation_cycle_list::iterator cIter = p_d.begin(); 
					cIter != p_d.end(); ++cIter) {
				
				/* sign-doubled cycle */
				permutation_cycle c_d = *cIter;
				
				/* Ignore identity cycle */
				if ( c_d.size() == 0 ) continue;
				
				/* Ignore cycles that start with a positive value - their
				 * negative complement is somewhere else in the cycle list
				 * (unless they're of the form (x,-x), in which case we want to
				 * ignore them anyway */
				if ( (c_d[0] & 0x1) == 0x0 ) continue;

				/* cycle in the original space */
				permutation_cycle c;
				permlib::dom_int x = c_d[0] >> 1U;
				c.push_back(x);

				/* Ensure that all values in this cycle are signed the same */
				for (uind i = 1; i < c_d.size(); ++i) {
					if ( (c_d[i] & 0x1) == 0x0 ) {
						valid = false;
						break;
					}

					/* removing the low order bit takes a row to it's original
					 * representation in this sign-doubling scheme */
					x = c_d[i] >> 1U;
					c.push_back(x);
				}

				/* stop looking at cycles for invalid permutations */
				if ( ! valid ) break;

				/* add un-doubled cycle to un-doubled permutation */
				p.push_back(c);
			}
			
			/* Skip invalid sign-doubled permutations */
			if ( ! valid ) continue;

//std::cout << **iter << "\n=> " << perm(g.dim(), p) << std::endl;

			/* add permutation to generator list */
			gens.push_back(boost::make_shared<permutation>(perm(g.dim(), p)));

		}

		return permlib::construct(g.dim(), gens.begin(), gens.end());
	}
	
	permutation_group_ptr shrink_group_to(permutation_group const& g, uind n) {
		/* generators of the new permutation group */
		std::vector<permutation_ptr> gens;

		permutation_list gen_list = strong_gen_set(g);
		for (permutation_list::iterator iter = gen_list.begin();
				iter != gen_list.end(); ++iter) {

			permutation_cycle_list p = cycle_list(**iter);
			permutation_cycle_list pn;

			/* valid permutations do not map i to j, j <= n < i. Since the
			 * cycles are sorted by starting element, it suffices to check that
			 * cycles starting with j <= n do not contain any elements i > n */
			for (permutation_cycle_list::iterator cIter = p.begin();
					cIter != p.end(); ++cIter) {

				permutation_cycle c = *cIter;

				/* Ignore identity cycle */
				if ( c.size() == 0 ) continue;

				/* stop when the cycles start to be > n (adjusted for 0-indexed
				 * cycles of 1-indexed permutations ...) */
				if ( c[0] >= n ) break;

				/* Ensure that all values in this cycle are <= n */
				bool valid = true;
				for (uind i = 1; i < c.size(); ++i) {
					if ( c[i] >= n ) {
						valid = false;
						break;
					}
				}
				if ( valid ) { pn.push_back(c); }
			}

			if ( pn.size() > 0 ) {
				gens.push_back(boost::make_shared<permutation>(perm(n, pn)));
			}
		}

		return permlib::construct(n, gens.begin(), gens.end());
	}

} /* namespace basil */
