/***********************************************************************
DESCRIPTION:
   Byte swapping routines used in various plugins
   There are two versions of each routine, one that's safe to use in
   all cases (but is slow) and one that is only safe to use on memory
   addresses that are aligned to the word size that's being byte-swapped
   but are much much much faster. Use the aligned versions of these
   routines whenever possible. The 'ndata' length count parameters and
   internal loops should be safe to use on huge memory arrays on 64-bit
   machines.

Notes:
    - adopted from VMD source code with some modifications

Revisions:
2012/03/13:
    - use system-independent types in the definitions of all functions
            size_t        <-    long     (for array index)
            int8_t        <-    char
            int16_t       <-    short
            int32_t       <-    int
***********************************************************************/

#ifndef __GEM_LIB_EDIAN_HPP__
#define __GEM_LIB_EDIAN_HPP__

#include <cstring>
#include <stdint.h>

namespace gem {

/* works on unaligned 2-byte quantities */
void swap2_unaligned(void *v, size_t ndata);

/* works on unaligned 4-byte quantities */
void swap4_unaligned(void *v, size_t ndata);

/* works on unaligned 8-byte quantities */
void swap8_unaligned(void *v, size_t ndata);

/* Only works with aligned 2-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap2_aligned(void *v, size_t ndata);

/* Only works with aligned 4-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap4_aligned(void *v, size_t ndata);

/* Only works with aligned 8-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap8_aligned(void *v, size_t ndata);

} // namespace gem

#endif
