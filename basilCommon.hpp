#ifndef _BASIL_COMMON_HPP_
#define _BASIL_COMMON_HPP_

#include <boost/shared_ptr.hpp>

#include <permlib/common.h> //because the PermLib author didn't ...
#include <permlib/bsgs.h>
#include <permlib/permutation.h>
#include <permlib/transversal/schreier_tree_transversal.h>

#include "lrs/lrs.hpp"
#include "lrs/matrix.hpp"

/** Namespace for the basil project.
 *  Imports types from boost and PermLib, and should have no naming conflicts 
 *  with the std namespace.
 */
namespace basil {
	/** import boost shared pointer into this namespace */
	using boost::shared_ptr;
	

	/** matrix type */
	typedef 
		lrs::matrix
		matrix;
	typedef
		shared_ptr<matrix>
		matrix_ptr;
	
	/** typesafe index into matrix */
	typedef 
		lrs::ind
		ind;
	/** unsigned version of ind */
	typedef
		lrs::uind
		uind;
	
	/** permutation type */
	typedef 
		permlib::Permutation 
		permutation;
	typedef
		shared_ptr<permutation>
		permutation_ptr;
		
	/** permutation tree traversal type */
	typedef 
		permlib::SchreierTreeTransversal<permutation>
		permutation_transversal;
	typedef
		shared_ptr<permutation_transversal>
		permutation_transversal_ptr;
	
	/** permutation group type */
	typedef 
		permlib::BSGS<permutation, permutation_transversal> 
		permutation_group;
	typedef
		shared_ptr<permutation_group>
		permutation_group_ptr;
	
	/** list of permutation type */
	typedef
		typename permutation_group::PERMlist
		permutation_list;
	
} /* namespace basil */

#endif /* _BASIL_COMMON_HPP_ */