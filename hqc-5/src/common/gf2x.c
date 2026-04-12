#include "gf2x.h"
#include "parameters.h"
#include <stddef.h>
#include <stdint.h>
#include "reduce.h"
#include "gfmul_fft.h"
#include <string.h>

/**
 * @brief Multiply two polynomials modulo \f$ X^n - 1\f$.
 *
 * This functions multiplies polynomials <b>v1</b> and <b>v2</b>.
 * The multiplication is done modulo \f$ X^n - 1\f$.
 *
 * @param[out] o Product of <b>v1</b> and <b>v2</b>
 * @param[in] v1 Pointer to the first polynomial
 * @param[in] v2 Pointer to the second polynomial
 */
void PQCLEAN_HQC256_CLEAN_vect_mul(uint64_t *o, const uint64_t *v1, const uint64_t *v2) {

	uint32_t o_u32[FFT_N / 32] = { 0 };

	bmul2_8192_to_16384((uint8_t*)o_u32, (uint8_t*)v1, (uint8_t*)v2);
	reduce((uint32_t*)o, (uint32_t*)o_u32);

}
void PQCLEAN_HQC256_CLEAN_vect_mul_2(uint64_t *u, uint64_t* tmp2, const uint64_t *r2, const uint64_t *h, const uint64_t *s)
{
	uint32_t h_u32[131072 / 32];
	uint32_t s_u32[131072 / 32];
	uint32_t r2_u32[131072 / 32];

	bmul2_8192_to_16384_prepare(h_u32, (uint8_t*)h);
	bmul2_8192_to_16384_prepare(r2_u32, (uint8_t*)r2);
	bmul2_8192_to_16384_prepare(s_u32, (uint8_t*)s);
	uint32_t o_u32[FFT_N / 32] = { 0 };

	bmul2_8192_to_16384_mul((uint8_t*)o_u32, r2_u32, h_u32);
	reduce((uint32_t*)u, (uint32_t*)o_u32);


	bmul2_8192_to_16384_mul((uint8_t*)o_u32, r2_u32, s_u32);
	reduce((uint32_t*)tmp2, (uint32_t*)o_u32);
}
