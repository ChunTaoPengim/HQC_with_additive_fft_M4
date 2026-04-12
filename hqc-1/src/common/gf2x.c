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
void PQCLEAN_HQC128_CLEAN_vect_mul(uint64_t *o, const uint64_t *v1, const uint64_t *v2) {

	uint32_t o_u32[FFT_N / 32] = { 0 };

	crt_full((uint8_t*)o_u32, (const uint8_t*)v1, (const uint8_t*)v2);
	reduce((uint32_t*)o, (uint32_t*)o_u32);

}
void PQCLEAN_HQC128_CLEAN_vect_mul_2(uint64_t *u, uint64_t* tmp2, const uint64_t *r2, const uint64_t *h, const uint64_t *s){

	uint32_t h_fft[32768 / 32] = {0};
    uint32_t s_fft[32768 / 32] = {0};
    uint32_t r2_fft[32768 / 32] = {0};

	fafft_input(h_fft, (uint8_t*)h);
    fafft_input(s_fft, (uint8_t*)s);
    fafft_input(r2_fft, (uint8_t*)r2);
	uint32_t o_u32[FFT_N / 32] = { 0 };

	crt_combine((uint8_t*)o_u32, (uint8_t*)r2, (uint8_t*)h, r2_fft, h_fft);
	reduce((uint32_t*)u, (uint32_t*)o_u32);

	crt_combine((uint8_t*)o_u32, (uint8_t*)r2, (uint8_t*)s, r2_fft, s_fft);
	reduce((uint32_t*)tmp2, (uint32_t*)o_u32);

}