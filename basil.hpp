#ifndef _BASIL_HPP_
#define _BASIL_HPP_

/** Namespace and common imports for the Basil project.
 *
 *  @author Aaron Moss
 */

/*  Copyright: Aaron Moss, 2012, moss.aaron@unb.ca  */

/*  This file is part of Basil.

    Basil is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Basil is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Basil.  If not, see <http://www.gnu.org/licenses/>.  */

#include <string>

#include <boost/shared_ptr.hpp>

#include "lrs/cobasis.hpp"
#include "lrs/matrix.hpp"

#include "permlib/common.h" //because the PermLib author didn't ...
#include "permlib/bsgs.h"
#include "permlib/permutation.h"
#include "permlib/transversal/schreier_tree_transversal.h"

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
	
	/** ordered list of matrix indices */
	typedef std::vector<permlib::dom_int> index_list;

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
