#ifndef _BASIL_COMMON_HPP_
#define _BASIL_COMMON_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gmpxx.h>

#include <permlib/common.h> //because the PermLib author didn't ...
#include <permlib/bsgs.h>
#include <permlib/permutation.h>
#include <permlib/transversal/schreier_tree_transversal.h>

namespace basil {
	//import boost shared pointer into this namespace
	using boost::shared_ptr;
	
	//value type of matrix
	typedef 
		mpz_class
		val_type;
	
	//matrix type
	typedef 
		boost::numeric::ublas::matrix<val_type> 
		matrix;
	//typesafe index into matrix
	typedef 
		matrix::size_type 
		ind;
	
	//permutation type
	typedef 
		permlib::Permutation 
		permutation;
	//permutation tree traversal type
	typedef 
		permlib::SchreierTreeTransversal<permutation>
		permutation_transversal;
	//permutation group type
	typedef 
		permlib::BSGS<permutation, permutation_transversal> 
		permutation_group;
}

#endif /* _BASIL_COMMON_HPP_ */