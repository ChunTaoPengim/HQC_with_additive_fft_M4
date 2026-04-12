#include "gf2x.h"
#include "parameters.h"
#include <stddef.h>
#include <stdint.h>

#define REDUCTION_MASKING (1<<(N&0x1F))-1

void reduce(uint32_t *o, const uint32_t *a) {
	uint32_t r;
	uint32_t carry;

	for (size_t i = 0; i < VEC_N_SIZE_32; ++i) {
		r = a[i + VEC_N_SIZE_32 - 1] >> (N & 0x1F);
		carry = a[i + VEC_N_SIZE_32] << (32 - (N & 0x1F));
		o[i] = a[i] ^ r ^ carry;
	}

	o[VEC_N_SIZE_32 - 1] &= REDUCTION_MASKING;
}