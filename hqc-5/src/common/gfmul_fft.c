// Implemented by Ming-Shing Chen, Tung Chou and Markus Krausz.
// 
// Modification :  Myeonghoon Lee and Jihoon Jang
//
// public domain

#include "parameters.h"
#include <stddef.h>
#include <stdint.h>
#include "encode.h"
#include "btfy_ffft.h"
#include "bitconversion.h"
#include "gfmul_fft.h"

#include <string.h>

//new code block
void bmul2_8192_to_16384_prepare(uint32_t * a_out, const uint8_t * a_in)
{
	uint32_t a0[65536 / 32] = { 0 };
	memcpy(a0, a_in, VEC_N_SIZE_64 * 8);
	bc_1_65536(a0);



	uint32_t temp[16];
	for (int i = 0; i < ((131072 / 32) / 32); i++) {
		for (int j = 0; j < 16; j++) temp[j] = a0[((131072 / 32) / 32)*j + i];
		encode_trunc_fft(a_out + i * 32, temp);
	}
	btfy_131072(a_out);
}


//new code block
void bmul2_8192_to_16384_mul(uint8_t * c, const uint32_t * a, const uint32_t * b)
{
	uint32_t a0[131072 / 32];
	uint32_t * c32 = (uint32_t*)c;

	for (int i = 0; i < (131072 / 32); i += 32) { gf232v_mul(a0 + i, a + i, b + i); }

	ibtfy_131072(a0);

	uint32_t temp[32];
	for (int i = 0; i < (131072 / 32); i += 32) {
		decode_trunc_fft(temp, a0 + i);
		for (int j = 0; j < 32; j++) c32[((131072 / 32) / 32)*j + i / 32] = temp[j];
	}


	ibc_1_131072(c32);

}

//new code block
void bmul2_8192_to_16384(uint8_t * c, const uint8_t * a, const uint8_t * b)
{
	uint32_t a_u32[131072 / 32];
	uint32_t b_u32[131072 / 32];

	bmul2_8192_to_16384_prepare(a_u32, a);
	bmul2_8192_to_16384_prepare(b_u32, b);

	bmul2_8192_to_16384_mul(c, a_u32, b_u32);
}
