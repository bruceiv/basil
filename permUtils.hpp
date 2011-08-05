#ifndef _PERM_UTILS_HPP_
#define _PERM_UTILS_HPP_

#include <set>
#include <sstream>

#include <permlib/permlib_api.h>

#include "basil.hpp"

namespace basil {
	
	/** Converts a permutation to a string in a cycle representation compatible 
	 *  with its input format. Derived from standard printing code in PermLib
	 *  @param p		The perumation to print
	 *  @return A representation of p as comma-delimited cycles of whitespace-
	 *  		delimited elements, where the cycles are ordered by minimum 
	 *  		element
	 */
	string in_str(permutation const& p) {
		typedef permlib::dom_int dom_int;
		std::ostringstream o;
		
		/* set of group elements that have been added to their cycle */
		std::set<dom_int> worked;
		/* for each element in the permutation */
		for (dom_int x = 0; x < p.size(); ++x) {
			dom_int px, startX;
			startX = x;		/* the start of the cycle */
			px = p/x;		/* the permutation applied to x */
			/* if x has already been found, or is not moved, ignore */
			if (worked.count(x) || x == px) {
				continue;
			}
			
			/* start the cycle */
			if ( ! worked.empty() ) o << " ,";
			worked.insert(x);
			o << " " << (x+1);
			/* while the cycle is not complete, follow it */
			while (p/px != startX) {
				o << " " << (px+1);
				worked.insert(px);
				px = p/px;
			}
			worked.insert(px);
			o << " " << (px+1);
		}
		return o.str();
	}
	
	/** Gets the strong generating set for a permutation group.
	 *  @param g		The permutation group
	 *  @return The strong generating set of generators for the group, as a 
	 *  		list of pointers to permutations.
	 */
	permutation_list strong_gen_set(permutation_group const& g) {
		return g.S;
	}
	
	/** Gets the strong generating set for a permutation group.
	 *  @param g		The permutation group
	 *  @return The strong generating set of generators for the group, as a 
	 *  		list of pointers to permutations.
	 */
	permutation_list& strong_gen_set(permutation_group& g) {
		return g.S;
	}
	
	/** Gets a small set of generators for a permutation group.
	 *  @param g		The permutation group
	 *  @return A hopefully small list of generators for the group, as a list 
	 *  		of pointers to permutations.
	 */
	permutation_list small_gen_set(permutation_group& g) {
		typedef permutation_list::iterator perm_iter;
		
		perm_iter iter = g.S.begin();
		permutation_list gens;
		if ( iter == g.S.end() ) return gens;
		
		/* target order */
		uint64_t ord = g.order();
		
		/* first, sort through and find out which generators are essential */
		permutation_list gsgs(g.S);
		permutation_list opts;
		
		iter = gsgs.begin();
		perm_iter next_iter = iter; ++next_iter;
		permutation_ptr p;
		permutation_group_ptr gn;
		
		while ( iter != gsgs.end() ) {
			/* remove a generator from the generator list */
			p = *iter;
			gsgs.erase(iter);
			
			/* construct a new group without it */
			gn = permlib::construct(g.n, gsgs.begin(), gsgs.end());
			
			/* test group order, keeping any permutation that lowers it */
			if ( gn->order() < ord ) gens.push_back(p); else opts.push_back(p);
			/* put permutation back in list regardless */
			gsgs.insert(next_iter, p);
			
			iter = next_iter; ++next_iter;
		}
		
		/* new group generated iteratively from the generators of g, initially 
		 * constructed to include all the essential generators */
		gn = permlib::construct(g.n, gens.begin(), gens.end());
		iter = opts.begin();
		
		/* while there are more possible generators, and the new group is not 
		 * of full order, add more generators */
		while ( iter != opts.end() && gn->order() < ord ) {
			/* if the current generator is not a member of the new group, 
			 * add it to the generators */
			if ( ! gn->sifts(**iter) ) {
				gens.push_back(*iter);
				gn = permlib::construct(g.n, gens.begin(), gens.end());
			}
			
			++iter;
		}
		
		return gens;
	}
	
} /* namespace basil */

#endif /* _PERM_UTILS_HPP_ */
