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
#include "bc_1.h"
#include "gfmul_fft.h"
#include "radix_mul.h"

#include <string.h>

#include "gen/IL_defs.h"
#include "gen/minpoly_mods_defs.h"

void fafft_input(uint32_t * a_out, const uint8_t * a_in)
{
	uint32_t a0[65536 / 32] = { 0 };
	memcpy(a0, a_in, VEC_N_SIZE_BYTES);
	bc_1_32768_plus_4096(a0);

	uint32_t temp[32];
	for (int i = 0; i < ((65536 / 32) / 32); i++) {
		for (int j = 0; j < 32; j++) temp[j] = a0[((65536 / 32) / 32)*j + i];
		encode_to_gft_full_length(a_out + i * 32, temp);
	}
	btfy_65536(a_out);
}

void fafft_output(uint8_t * c, const uint32_t * a, const uint32_t * b)
{
	uint32_t a0[65536 / 32];
	uint32_t * c32 = (uint32_t*)c;

	for (int i = 0; i < (65536 / 32); i += 32) { gf232v_mul(a0 + i, a + i, b + i); }

	ibtfy_65536(a0);

	uint32_t temp[32];
	for (int i = 0; i < (65536 / 32); i += 32) {
		decode_from_gft_full_length(temp, a0 + i);
		for (int j = 0; j < 32; j++) c32[((65536 / 32) / 32)*j + i / 32] = temp[j];
	}

	ibc_1_65536(c32);
}

void crt_combine(uint8_t * c, const uint8_t * a, const uint8_t * b, const uint32_t * a_fft, const uint32_t * b_fft)
{
	uint32_t pc_m65536[65536 / 32];
	uint32_t pc_m8192[8192 / 32];
	fafft_output((uint8_t*)pc_m65536, a_fft, b_fft);
	uint32_t temp2[256];
	gf2x_mul_8192_no_inv_swap(pc_m8192, (uint32_t*)a, (uint32_t*)b);

	poly_to_rd16(temp2, pc_m65536, 256);
	// the mod x^8192
	for (int i = 0; i < (8192 / 32); i++) {
		pc_m8192[i] ^= temp2[i];
	}

	// build pc_mix * IL
	gf2x_mul_8192_preswap(pc_m8192, IL_swapped_65536, pc_m8192);

	uint32_t pc_icrt0[(65536+8192) / 32] __attribute__((aligned(64))) = { 0 };

	// sparse polynomial
    {
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < minpoly_mod_lens[0]; j++) {
                pc_icrt0[i + minpoly_mods[0][j]] ^= pc_m8192[i];
            }
        }
    }
    #pragma GCC unroll 31
    for (int rem = 1; rem < 32; rem++) {
        for (int i = 0; i <= 256; i++) {
            uint32_t value = (i > 0) ? (pc_m8192[i - 1] >> (32 - rem)) : 0;
            if (i < 256) value |= pc_m8192[i] << rem;

            for (int j = 0; j < minpoly_mod_lens[rem]; j++) {
                pc_icrt0[i + minpoly_mods[rem][j]] ^= value;
            }
        }
    }

	for(int i = 0; i < (65536)/32; i++)
	{
		pc_icrt0[i] ^= pc_m65536[i];
	}

	memcpy(c, pc_icrt0, sizeof(uint32_t)*((65536+8192) / 32));
}

void crt_full(uint8_t * c, const uint8_t * a, const uint8_t * b)
{

	uint32_t a_u32[65536 / 32];
	uint32_t b_u32[65536 / 32];
	uint32_t pc_m65536[65536 / 32];
	uint32_t pc_m8192[8192 / 32];
	fafft_input(a_u32, a);
	fafft_input(b_u32, b);
	fafft_output((uint8_t*)pc_m65536, a_u32, b_u32);
	uint32_t temp2[256];
	gf2x_mul_8192_no_inv_swap(pc_m8192, (uint32_t*)a, (uint32_t*)b);

	poly_to_rd16(temp2, pc_m65536, 256);
	// the mod x^8192
	for (int i = 0; i < (8192 / 32); i++) {
		pc_m8192[i] ^= temp2[i];
	}

	// build pc_mix * IL
	gf2x_mul_8192_preswap(pc_m8192, IL_swapped_65536, pc_m8192);
	uint32_t pc_icrt0[(65536+8192) / 32] __attribute__((aligned(64))) = { 0 };

	// sparse polynomial
    {
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < minpoly_mod_lens[0]; j++) {
                pc_icrt0[i + minpoly_mods[0][j]] ^= pc_m8192[i];
            }
        }
    }
    #pragma GCC unroll 31
    for (int rem = 1; rem < 32; rem++) {
        for (int i = 0; i <= 256; i++) {
            uint32_t value = (i > 0) ? (pc_m8192[i - 1] >> (32 - rem)) : 0;
            if (i < 256) value |= pc_m8192[i] << rem;

            for (int j = 0; j < minpoly_mod_lens[rem]; j++) {
                pc_icrt0[i + minpoly_mods[rem][j]] ^= value;
            }
        }
    }

	for(int i = 0; i < (65536)/32; i++)
	{
		pc_icrt0[i] ^= pc_m65536[i];
	}

	memcpy(c, pc_icrt0, sizeof(uint32_t)*((65536+8192) / 32));
}
