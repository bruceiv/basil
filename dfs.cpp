#include "basilCommon.hpp"
#include "dfs.hpp"
#include "lrs/lrs.hpp"

namespace basil {
	
	void dfs::doDfs() {
		if (opts.showAllDicts) l.printDict();
		
		l.getFirstBasis();
		
		if (opts.showAllDicts) l.printDict();
		
	}

	
} /* namespace basil */
