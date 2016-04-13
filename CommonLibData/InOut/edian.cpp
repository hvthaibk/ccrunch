#include "edian.hpp"

namespace gem {

/* works on unaligned 2-byte quantities */
void swap2_unaligned(void *v, size_t ndata)
{
    size_t    i;
    int8_t    *dataptr = (int8_t *) v;
    int8_t    tmp;

    for (i = 0; i < ndata-1; i += 2) {
        tmp = dataptr[i];
        dataptr[i] = dataptr[i+1];
        dataptr[i+1] = tmp;
    }
}

/* works on unaligned 4-byte quantities */
void swap4_unaligned(void *v, size_t ndata)
{
    size_t    i;
    int8_t    *dataptr = (int8_t *) v;
    int8_t    tmp;

    for (i = 0; i < ndata; i++) {
        tmp = dataptr[0];
        dataptr[0] = dataptr[3];
        dataptr[3] = tmp;
        tmp = dataptr[1];
        dataptr[1] = dataptr[2];
        dataptr[2] = tmp;
        dataptr += 4;
    }
}

/* works on unaligned 8-byte quantities */
void swap8_unaligned(void *v, size_t ndata)
{
    size_t    i;
    int8_t    *data = (int8_t *) v;
    int8_t    byteArray[8];
    int8_t    *bytePointer;

    for (i = 0; i < ndata; i++) {
        bytePointer = data + (i<<3);
        byteArray[0]  =  *bytePointer;
        byteArray[1]  =  *(bytePointer+1);
        byteArray[2]  =  *(bytePointer+2);
        byteArray[3]  =  *(bytePointer+3);
        byteArray[4]  =  *(bytePointer+4);
        byteArray[5]  =  *(bytePointer+5);
        byteArray[6]  =  *(bytePointer+6);
        byteArray[7]  =  *(bytePointer+7);

        *bytePointer     = byteArray[7];
        *(bytePointer+1) = byteArray[6];
        *(bytePointer+2) = byteArray[5];
        *(bytePointer+3) = byteArray[4];
        *(bytePointer+4) = byteArray[3];
        *(bytePointer+5) = byteArray[2];
        *(bytePointer+6) = byteArray[1];
        *(bytePointer+7) = byteArray[0];
    }
}

/* Only works with aligned 2-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap2_aligned(void *v, size_t ndata)
{
    size_t     i;
    int16_t    *data = (int16_t *) v;
    int16_t    *N;

    for (i = 0; i < ndata; i++) {
        N  = data + i;
        *N = (int16_t) (((*N>>8)&0xff) | ((*N&0xff)<<8));
    }
}

/* Only works with aligned 4-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap4_aligned(void *v, size_t ndata)
{
    size_t     i;
    int32_t    *data = (int32_t *) v;
    int32_t    *N;

    for (i = 0; i < ndata; i++) {
        N  = data + i;
        *N = (((*N>>24)&0xff) | ((*N&0xff)<<24) |
             ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
    }
}

/* Only works with aligned 8-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
void swap8_aligned(void *v, size_t ndata)
{
    /* Use int32_t* internally to prevent bugs caused by some compilers */
    /* and hardware that would potentially load data into an FP reg     */
    /* and hose everything, such as the old "jmemcpy()" bug in NAMD     */
    size_t     i;
    int32_t    *data = (int32_t *) v;
    int32_t    *N;
    int32_t    t0, t1;

    for (i = 0; i < ndata; i++) {
        N  = data + (i<<1);
        t0 = N[0];
        t0 = (((t0>>24)&0xff) | ((t0&0xff)<<24) |
             ((t0>>8)&0xff00) | ((t0&0xff00)<<8));

        t1 = N[1];
        t1 = (((t1>>24)&0xff) | ((t1&0xff)<<24) |
             ((t1>>8)&0xff00) | ((t1&0xff00)<<8));

        N[0] = t1;
        N[1] = t0;
    }
}

} // namespace gem
