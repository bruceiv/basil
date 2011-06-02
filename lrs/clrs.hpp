#ifndef _CLRS_HPP_
#define _CLRS_HPP_

#include "lrslib.h"
#include "lrsgmp.h"
// because the LRS author saw fit to polute the global namespace with this ....
#undef copy

lrs_mp_matrix lrs_alloc_mp_matrix(long,long);
void lrs_clear_mp_matrix(lrs_mp_matrix,long,long);

#endif /* _CLRS_HPP_ */