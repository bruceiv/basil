#ifndef _DFS_HPP_
#define _DFS_HPP_

#include "basilCommon.hpp"
#include "lrs/lrs.hpp"

namespace basil {
	
	/** Stateful wrapper class for DFS algorithm.
	 */
	class dfs {
	public:
		/** DFS on the given matrix, according to the provided permutation group.
		 *  @param m		The matrix to DFS on
		 *  @param g		The permutation group of the matrix
		 */
		void doDfs(matrix m, permutation_group g);
	
	protected:
		/** LRS wrapper for this DFS */
		lrs::lrs l;
	}; /* class dfs */

} /* namespace basil */
#endif /* _DFS_HPP_ */