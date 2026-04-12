#ifndef _GFMUL_FFT_H_
#define _GFMUL_FFT_H_

#include <stdint.h>

void fafft_input(uint32_t * a_out, const uint8_t * a_in);
void fafft_output(uint8_t * c, const uint32_t * a, const uint32_t * b);
void crt_combine(uint8_t * c, const uint8_t * a, const uint8_t * b, const uint32_t * a_fft, const uint32_t * b_fft);
void crt_full(uint8_t * c, const uint8_t * a, const uint8_t * b);

#endif
