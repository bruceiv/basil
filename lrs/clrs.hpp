#ifndef _CLRS_HPP_
#define _CLRS_HPP_

extern "C" {
#include "lrslib.h"
#include "lrsgmp.h"
// because the LRS author saw fit to polute the global namespace with this ....
#undef copy
}

#endif /* _CLRS_HPP_ */