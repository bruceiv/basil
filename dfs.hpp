#ifndef _DFS_HPP_
#define _DFS_HPP_

#include "basilCommon.hpp"
#include "lrs/lrs.hpp"

namespace basil {
	
	/** Options for DFS algorithm.
	 *  The default constructor initializes the data members to their default 
	 *  values, while the methods can be used to modify the options.
	 */
	struct dfs_opts {
		
		/** Default constructor - sets all options to their default value */
		dfs_opts() : showAllDicts(false) {}
		
		/** Activates (or deactivates) the showsAllDicts option */
		dfs_opts& showsAllDicts(bool opt = true) { showAllDicts = opt;}
		
		/** show all dictionaries as they are generated */
		bool showAllDicts;
	}; /* struct dfs_opts */
	
	/** Stateful wrapper class for DFS algorithm. */
	class dfs {
	public:
		/** Set up a DFS on the given matrix, with the given permuation group.
		 *  @param m		The matrix to DFS on
		 *  @param g		The permutation group of the matrix
		 *  @param opts		The options for this DFS (default values if not 
		 * 					provided)
		 */
		dfs(matrix& m, permutation_group& g, dfs_opts opts = dfs_opts()) 
				: l(m), g(g), opts(opts) {}
		
		/** Perform the DFS. */
		void doDfs();
	
	protected:
		/** LRS wrapper for this DFS */
		lrs::lrs l;
		/** Permutation group used for this DFS */
		permutation_group& g;
		/** Options for controlling the DFS algorithm */
		dfs_opts opts;
	}; /* class dfs */

} /* namespace basil */
#endif /* _DFS_HPP_ */