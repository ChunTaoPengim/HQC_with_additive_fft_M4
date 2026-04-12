#ifndef _RADIX_MUL_H
#define _RADIX_MUL_H

#include <stdint.h>

void gf2x_mul_3072(uint32_t *c, const uint32_t *a, const uint32_t *b);
void gf2x_mul_3072_no_inv_swap(uint32_t *c, const uint32_t *a, const uint32_t *b);
void gf2x_mul_3072_preswap(uint32_t *c, const uint32_t *a, const uint32_t *b);
void poly_to_rd16(uint32_t *res, const uint32_t * x, uint32_t len);
void rd16_to_poly(uint32_t * x, uint32_t len);

#endif
