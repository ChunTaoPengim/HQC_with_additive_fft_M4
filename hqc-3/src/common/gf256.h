#ifndef _GF256_H
#define _GF256_H

#include <stdint.h>

void gf256_madd(uint32_t *c, const uint32_t *a, const uint32_t b);
void gf256_mul(uint32_t *a, uint32_t b);
void gf256_madd_16(uint32_t *c, const uint32_t *a, uint32_t b);
uint32_t radix_16_gf256_4x4(uint32_t a, uint32_t b);
uint32_t radix_16_matrix_vector(uint32_t *mat, uint32_t *vec, int vec_len);
uint32_t radix_16_matrix_vector_half(uint32_t *mat, uint32_t *vec, int vec_len);

#endif
