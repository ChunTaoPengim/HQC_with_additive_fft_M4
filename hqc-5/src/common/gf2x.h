#ifndef GF2X_H
#define GF2X_H
/**
 * @file gf2x.h
 * @brief Header file for gf2x.c
 */

#include <stdint.h>

void PQCLEAN_HQC256_CLEAN_vect_mul(uint64_t *o, const uint64_t *v1, const uint64_t *v2);
void PQCLEAN_HQC256_CLEAN_vect_mul_2(uint64_t *u, uint64_t* tmp2, const uint64_t *r2, const uint64_t *h, const uint64_t *s);

#endif
