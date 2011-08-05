#include <boost/make_shared.hpp>

#include <permlib/common.h>
#include <permlib/sorter/base_sorter.h>
//above included because the PermLib author didn't ...
#include <permlib/symmetric_group.h>
#include <permlib/search/partition/matrix_automorphism_search.h>

#include "automorphism.hpp"
#include "basil.hpp"
#include "gram.hpp"

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
		
		/* the search to perform */
		matrix_automorphism_search mas(s_n, 0);
		mas.construct(g);
		
		/* perform search */
		permutation_group_ptr ret = boost::make_shared<permutation_group>(n);
		mas.search(*ret);
		
		/* return results */
		return ret;
	}
	
} /* namespace basil */
