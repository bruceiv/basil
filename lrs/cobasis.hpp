#ifndef _COBASIS_H_
#define _COBASIS_H_

#include <vector>

#include <gmpxx.h>


namespace lrs {
	
	/** Represents a list of indexes */
	typedef std::vector<long> index_list;
	
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
