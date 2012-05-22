#ifndef _BASIL_HPP_
#define _BASIL_HPP_

#include <string>

#include <boost/shared_ptr.hpp>

#include <permlib/common.h> //because the PermLib author didn't ...
#include <permlib/bsgs.h>
#include <permlib/permutation.h>
#include <permlib/transversal/schreier_tree_transversal.h>

#include "lrs/cobasis.hpp"
#include "lrs/matrix.hpp"

namespace basil {
	
	////////////////////////////////////////////////////////////////////////////
	//
	//  Imports and typedefs for use in Basil
	//
	////////////////////////////////////////////////////////////////////////////
	
	/** import STL string into this namespace */
	using std::string;
	
	/** import boost shared pointer into this namespace */
	using boost::shared_ptr;
	

	/** matrix type */
	typedef 
		lrs::matrix_mpq
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
	
	/** efficient set of matrix indices */
	typedef 
		lrs::index_set 
		index_set;
	typedef 
		shared_ptr<index_set> 
		index_set_ptr;
	
	/** permutation type */
	typedef 
		permlib::Permutation 
		permutation;
	typedef
		shared_ptr<permutation>
		permutation_ptr;
		
	/** permutation traversal type */
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
		permutation_group::PERMlist
		permutation_list;

} /* namespace basil */

#endif /* _BASIL_HPP_ */
