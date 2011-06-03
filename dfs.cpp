#include "basilCommon.hpp"
#include "dfs.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/lrs.hpp"

namespace basil {
	
	void dfs::doDfs() {
		using lrs::cobasis;
		
		if (opts.showAllDicts) l.printDict();
		
		l.getFirstBasis();
		
		if (opts.showAllDicts) l.printDict();
		
		cobasis* cob = l.getCobasis(0);
	}

	
} /* namespace basil */
