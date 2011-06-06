#ifndef _COBASIS_H_
#define _COBASIS_H_

#include <vector>

#include <gmpxx.h>

#include "clrs.hpp"

namespace lrs {
	
	struct cobasis {
		
		cobasis(mpz_class& det, long ray, index_list& cob, 
				long totalInc, index_list& extraInc)
				: det(det), ray(ray), cob(cob), totalInc(totalInc), 
				extraInc(extraInc) {}
		
		mpz_class det;
		long ray;
		index_list cob;
		long totalInc;
		index_list extraInc;
	};
	
}

#endif /* _COBASIS_H_ */
