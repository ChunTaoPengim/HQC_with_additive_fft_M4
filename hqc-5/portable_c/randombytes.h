#ifndef _RANDOMBYTES_H_
#define _RANDOMBYTES_H_

//#define randombytes     PQCLEAN_randombytes

#include "stdlib.h"

static inline
int randombytes(void *output, size_t n)
{
  unsigned char * ptr = output;
  for(unsigned i=0;i<n;i++) { ptr[i]=rand(); }
  return n;
}




#endif



