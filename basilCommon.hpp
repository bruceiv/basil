#ifndef _BASIL_COMMON_HPP_
#define _BASIL_COMMON_HPP_

#include <boost/shared_ptr.hpp>

#include <permlib/common.h> //because the PermLib author didn't ...
#include <permlib/bsgs.h>
#include <permlib/permutation.h>
#include <permlib/transversal/schreier_tree_transversal.h>

#include "lrs/lrs.hpp"

/** Namespace for the basil project.
 *  Imports types from boost and PermLib, and should have no naming conflicts 
 *  with the std namespace.
 */
namespace basil {
	/** import boost shared pointer into this namespace */
	using boost::shared_ptr;
	
	/** value type of matrix */
	typedef 
		lrs::val_t
		val_type;
	
	/** matrix type */
	typedef 
		lrs::matrix
		matrix;
	/** typesafe index into matrix */
	typedef 
		lrs::ind
		ind;
	
	/** permutation type */
	typedef 
		permlib::Permutation 
		permutation;
	/** permutation tree traversal type */
	typedef 
		permlib::SchreierTreeTransversal<permutation>
		permutation_transversal;
	/** permutation group type */
	typedef 
		permlib::BSGS<permutation, permutation_transversal> 
		permutation_group;
} /* namespace basil */

#endif /* _BASIL_COMMON_HPP_ */