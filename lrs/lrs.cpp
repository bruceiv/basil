/* implements lrs.hpp C++ wrapper for LRS */

#include "lrs.hpp"

#include "lrslib.h"


namespace lrs {
	
	lrs::lrs() {
		lrs_init((char*)"basil");
	}
	
	lrs::~lrs() {
		lrs_close((char*)"basil");
	}
	
} /* namespace lrs */
