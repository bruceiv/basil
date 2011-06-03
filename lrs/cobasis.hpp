#ifndef _COBASIS_H_
#define _COBASIS_H_

#include <valarray>

#include <gmpxx.h>

namespace lrs {
	
	class cobasis {
	public:
		
		cobasis(mpz_class det, long ray, std::valarray<long> cob, 
				long totalInc, std::valarray<long> extraInc)
				: det_(det), ray_(ray), cob_(cob), totalInc_(totalInc), 
				extraInc_(extraInc) {}
		
		const mpz_class& det() const { return det_; }
		const long& ray() const { return ray_; }
		const std::valarray<long>& cob() const { return cob_; }
		const long& totalInc() const { return totalInc_; }
		const std::valarray<long>& extraInc() const { return extraInc_; }
		
	private:
		mpz_class det_;
		long ray_;
		std::valarray<long> cob_;
		long totalInc_;
		std::valarray<long> extraInc_;
	};
	
}

#endif /* _COBASIS_H_ */
