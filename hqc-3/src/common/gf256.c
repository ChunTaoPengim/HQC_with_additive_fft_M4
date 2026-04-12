#include "gf256.h"
#include "parameters.h"
#include <stddef.h>
#include <stdint.h>
#include "gf256_rdx16.h"

#if !defined(_M4_ASM_)

static inline
void bs_gf256_multab( uint32_t multab[8] , uint8_t b )
{
    uint32_t bx = b;
    multab[0] = bx;
    for(int i=1;i<8;i++) {
        uint32_t t1 = (bx>>7)&1;
        t1 *= 0x11d;
        bx = (bx<<1) ^ t1;
        multab[i] = bx;
    }
}

static inline
uint32_t bs_gf256_mul( uint32_t a , const uint32_t multab[8] ) {
    uint32_t res = 0;
    for(int i=0;i<8;i++) {
        uint32_t ai = (a>>i)&0x01010101;
        res ^= multab[i]*ai;
    }
    return res;
}

/**
 * Implement GF256 multiply-accumulate on 36 elements packed into 9 registers.
 * Input: First multiplicand a (every 4 field elements packed into one register)
 * Input: Second multiplicand b (1 field element in the least significant byte)
 * Input/Output: accumulator c (every 4 field elements packed into one register)
 * 0x11d : 1_0001_1101 gf256
 * c <- c + a*b
 */
void gf256_madd(uint32_t *c, const uint32_t *a, uint32_t b) {
    uint32_t regs[9] = {0};
    uint32_t m0[8];
    bs_gf256_multab(m0, (uint8_t)b);
    for(int j=0;j<9;j++) { regs[j] ^= bs_gf256_mul( a[j], m0 ); }
    for(int k=0;k<9;k++) { c[k] ^= regs[k]; }
}

/**
 * Implement GF256 multiply-accumulate on 16 elements packed into 4 registers.
 * Input: First multiplicand a (every 4 field elements packed into one register)
 * Input: Second multiplicand b (1 field element in the least significant byte)
 * Input/Output: accumulator c (every 4 field elements packed into one register)
 * 0x11d : 1_0001_1101 gf256
 * c <- c + a*b
 */
void gf256_madd_16(uint32_t *c, const uint32_t *a, uint32_t b) {
    uint32_t regs[4] = {0};
    uint32_t m0[8];
    bs_gf256_multab(m0, (uint8_t)b);
    for(int j=0;j<4;j++) { regs[j] ^= bs_gf256_mul( a[j], m0 ); }
    for(int k=0;k<4;k++) { c[k] ^= regs[k]; }
}

/**
 * Implement GF256 multiplication on 36 elements packed into 9 registers.
 * Input/Output: First multiplicand a (every 4 field elements packed into one register)
 * Input: Second multiplicand b (1 field element in the least significant byte)
 * 0x11d : 1_0001_1101 gf256
 * a <- a*b
 */
void gf256_mul(uint32_t *a, uint32_t b) {
    uint32_t regs[9] = {0};
    uint32_t m0[8];
    bs_gf256_multab(m0, (uint8_t)b);
    for(int j=0;j<9;j++) { regs[j] ^= bs_gf256_mul( a[j], m0 ); }
    for(int k=0;k<9;k++) { a[k] = regs[k]; }
}

static inline uint64_t _load(uint32_t in0, uint32_t in1) { return in0 | (((uint64_t)in1) << 32); }
static inline void _store(uint32_t *out0, uint32_t *out1, uint64_t in) {  *out0 = (uint32_t)(in & 0xFFFFFFFF);  *out1 = (uint32_t)(in >> 32); }

uint32_t radix_16_gf256_4x4(uint32_t a, uint32_t b) {
    uint32_t k[4] = {0, 2, 1, 3};
    uint32_t t[4] = {0};
    uint32_t tmp1 = 0;
    uint32_t tmp2 = 0;
    for(int i=0; i<4; i++){
        uint64_t a_i = (a >> k[i]) & _NORM_RDX16_;
        uint64_t b_i = (b >> k[i]) & _NORM_RDX16_;
        _store(&t[i], &tmp1, a_i * b_i);
        tmp1 &= _NORM_RDX16_;
        _store(&tmp2, &tmp1, (uint64_t)tmp1 * (uint64_t)_0x1D_RDX16_);
        tmp2 += tmp1 * _0x1D_RDX16_;
        t[i] ^= tmp2;
        t[i] &= _NORM_RDX16_;
    }
    uint32_t out = 0;
    out = t[0] | (t[1] << 2);
    out |= (t[2] << 1);
    out |= (t[3] << 3);

    return out;
}

/**
 * Implement GF256 matrix-vector multiplication in radix-16 representation
 * Input: matrix: a pointer to a pre-computed (4 x vec_len) matrix in radix-16 representation
 * Input: vector: a vector of vec_len elements
 * Output: result: the product vector of length 4 in radix-16 representation
 * 0x11d : 1_0001_1101 gf256
 */
uint32_t radix_16_matrix_vector(uint32_t *mat, uint32_t *vec, int vec_len) {

    uint64_t t0 = 0;
    uint64_t t1 = 0;
    uint64_t t2 = 0;
    uint64_t t3 = 0;
    for(int j=0;j<vec_len;j++) {
        t0 = rdx16_muladd( t0 , (mat[j])&_NORM_RDX16_ , vec[j] ) & _U64NORM_RDX16_;
        t1 = rdx16_muladd( t1 , (mat[j]>>2)&_NORM_RDX16_ , vec[j] ) & _U64NORM_RDX16_;
        t2 = rdx16_muladd( t2 , (mat[j]>>1)&_NORM_RDX16_ , vec[j] ) & _U64NORM_RDX16_;
        t3 = rdx16_muladd( t3 , (mat[j]>>3)&_NORM_RDX16_ , vec[j] ) & _U64NORM_RDX16_;
    }        
    uint32_t r0 = (((uint32_t)t0) ^ rdx16_gf256_reduce( t0>>32 ))&_NORM_RDX16_;
    uint32_t r1 = (((uint32_t)t1) ^ rdx16_gf256_reduce( t1>>32 ))&_NORM_RDX16_;
    uint32_t r2 = (((uint32_t)t2) ^ rdx16_gf256_reduce( t2>>32 ))&_NORM_RDX16_;
    uint32_t r3 = (((uint32_t)t3) ^ rdx16_gf256_reduce( t3>>32 ))&_NORM_RDX16_;
    
    return ( r0 | (r1<<2) | (r2<<1) | (r3<<3) );
}

/**
 * Implement GF256 matrix-vector multiplication in radix-16 representation with half size
 * Input: matrix: a pointer to a pre-computed (4 × vec_len) matrix in radix-16 representation
 * Input: vector: a vector of vec_len elements
 * Output: result: the product vector of length 4 in radix-16 representation
 * 0x11d : 1_0001_1101 gf256
 */
uint32_t radix_16_matrix_vector_half(uint32_t *mat, uint32_t *vec, int vec_len) {

    uint64_t t0 = 0;
    uint64_t t1 = 0;
    for(int j=0;j<vec_len;j++) {
        t0 = rdx16_muladd( t0 , (mat[j])&_NORM_RDX16_ , vec[j] ) & _U64NORM_RDX16_;
        t1 = rdx16_muladd( t1 , (mat[j]>>2)&_NORM_RDX16_ , vec[j] ) & _U64NORM_RDX16_;
    }        
    uint32_t r0 = (((uint32_t)t0) ^ rdx16_gf256_reduce( t0>>32 ))&_NORM_RDX16_;
    uint32_t r1 = (((uint32_t)t1) ^ rdx16_gf256_reduce( t1>>32 ))&_NORM_RDX16_;
    
    return ( r0 | (r1<<2) );
}

#endif
