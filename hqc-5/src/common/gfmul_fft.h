#ifndef _GFMUL_FFT_H_
#define _GFMUL_FFT_H_

#include <stdint.h>

void bmul2_8192_to_16384_prepare(uint32_t * a_out, const uint8_t * a_in);
void bmul2_8192_to_16384_mul(uint8_t * c, const uint32_t * a, const uint32_t * b);
void bmul2_8192_to_16384(uint8_t * c, const uint8_t * a, const uint8_t * b);

#endif
