#include "parameters.h"
#include "radix_mul.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#if defined(_M4_ASM_)

#define SWAPMOVE_CONST(RES, A, mask, n, tmp) \
{ \
    asm volatile( \
    	"EOR.W %0, %1, %1, LSR #"#n"\n" \
	: "=r"(tmp) \
	: "r"(A));\
	tmp &= (mask); \
    RES = A ^ tmp;\
    asm volatile( \
    	"EOR.W %0, %0, %1, LSL #"#n"\n" \
	: "+r"(RES) \
	: "r"(tmp));\
}

#define SWAPMOVE(A, B, mask, n, tmp) \
{ \
    asm volatile( \
    	"EOR.W %0, %1, %2, LSR #"#n"\n" \
	: "=r"(tmp) \
	: "r"(B), "r"(A));\
	tmp &= (mask); \
    B ^= tmp; \
    asm volatile( \
    	"EOR.W %0, %0, %1, LSL #"#n"\n" \
	: "+r"(A) \
	: "r"(tmp));\
}

void radix_16_basemul_asm(uint32_t *out, uint32_t a, uint32_t b);
void mul2_asm(uint32_t *out, const uint32_t *a, const uint32_t *b);
uint32_t radix_16_truncate_asm(uint32_t a, uint32_t b);

#else

#define SWAPMOVE_CONST(RES, A, mask, n, tmp) do { \
    tmp = (A ^ (A >> n)) & (mask); \
    RES = A ^ tmp; \
    RES ^= (tmp << n); \
} while(0)

#define SWAPMOVE(A, B, mask, n, tmp) do { \
    tmp = (B ^ (A >> n)) & (mask); \
    B ^= tmp; \
    A ^= (tmp << n); \
} while(0)

static inline uint64_t _load(uint32_t in0, uint32_t in1) { return in0 | (((uint64_t)in1) << 32); }
static inline void _store(uint32_t *out0, uint32_t *out1, uint64_t in) {  *out0 = (uint32_t)(in & 0xFFFFFFFF);  *out1 = (uint32_t)(in >> 32); }

static inline void radix_16_basemul_asm(uint32_t *out, uint32_t a, uint32_t b)
{
    uint32_t mask32 = 0x11111111;
    uint64_t a0 = (a)&mask32;
    uint64_t b0 = (b)&mask32;
    uint64_t b1 = (b >> 2) & mask32;
    uint64_t b2 = (b >> 1) & mask32;
    uint64_t b3 = (b >> 3) & mask32;
    uint32_t tmp[8];
//uint32 0, [out[0], tmp[1]] = a0 * b0
    _store(&out[0], &tmp[1], a0 * b0);
    out[0] &= mask32;
//uint32 1, [tmp[2], tmp[3]] = a0 * b2
    _store(&tmp[2], &tmp[3], a0 * b2);
//uint32 2, [tmp[1], tmp[2]] += a0 * b1
    _store(&tmp[1], &tmp[2], _load(tmp[1], tmp[2]) + a0 * b1);
    tmp[1] &= mask32;
    uint64_t a1 = (a >> 2) & mask32;
//uint32_t 3, [tmp[1], tmp[2]] += a1 * b0;
    _store(&tmp[1], &tmp[2], _load(tmp[1], tmp[2]) + a1 * b0);
    tmp[1] &= mask32;
    tmp[2] &= mask32;
    out[0] ^= (tmp[1]<<2);
//uint32_t 4, [tmp[4], tmp[5]] = a1 * b3;
    _store(&tmp[4], &tmp[5], a1 * b3);
//uint32_t 5, [tmp[3], tmp[4]] += a0 * b3
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a0 * b3);
//uint32_t 6, [tmp[2], tmp[3]] += a1 * b1
    _store(&tmp[2], &tmp[3], _load(tmp[2], tmp[3]) + a1 * b1);
    tmp[2] &= mask32;
    tmp[3] &= mask32;
//uint32_t 7, [tmp[3], tmp[4]] += a1 * b2;
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a1 * b2);
    tmp[4] &= mask32;
    uint64_t a2 = (a >> 1) & mask32;
//uint32_t 8, [tmp[2], tmp[3]] += a2 * b0;
    _store(&tmp[2], &tmp[3], _load(tmp[2], tmp[3]) + a2 * b0);
    tmp[2] &= mask32;
    tmp[3] &= mask32;
    out[0] ^= (tmp[2]<<1);
//uint32_t 9, [tmp[3], tmp[4]] += a2 * b1
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a2 * b1);
    tmp[3] &= mask32;
//uint32_t 10, [tmp[4], tmp[5]] += a2 * b2;
    _store(&tmp[4], &tmp[5], _load(tmp[4], tmp[5]) + a2 * b2);
    tmp[4] &= mask32;
    uint64_t a3 = (a >> 3) & mask32;
//uint32_t 11, [tmp[3], tmp[4]] += a3 * b0;
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a3 * b0);
    tmp[3] &= mask32;
    out[0] ^= (tmp[3]<<3);
//uint32_t 12, [tmp[6], tmp[7]] = a3 * b3;
    _store(&tmp[6], &tmp[7], a3 * b3);
    tmp[7] &= mask32;
//uint32_t 13, [tmp[5], tmp[6]] += a2 * b3;
    _store(&tmp[5], &tmp[6], _load(tmp[5], tmp[6]) + a2 * b3);
    tmp[5] &= mask32;
//uint32_t 14, [tmp[4], tmp[5]] += a3 * b1;
    _store(&tmp[4], &tmp[5], _load(tmp[4], tmp[5]) + a3 * b1);
    tmp[4] &= mask32;
//uint32_t 15, [tmp[5], tmp[6]] += a3 * b2;
    _store(&tmp[5], &tmp[6], _load(tmp[5], tmp[6]) + a3 * b2);
    tmp[5] &= mask32;
    tmp[6] &= mask32;
    out[1] = tmp[4] ^ (tmp[5]<<2);
    out[1] ^= (tmp[6]<<1);
    out[1] ^= (tmp[7]<<3);
    out[0] = out[0];
    out[1] = out[1];
}

static void mul2_asm(uint32_t *out, const uint32_t *a, const uint32_t *b)
{
   uint32_t hs0, hs1;
   uint32_t hl2[2];

   hs0 = a[0] ^ a[1];
   hs1 = b[0] ^ b[1];

   radix_16_basemul_asm(out, a[0], b[0]);
   radix_16_basemul_asm(out+2, a[1], b[1]);
   radix_16_basemul_asm(hl2, hs0, hs1);


   hl2[0] = hl2[0] ^ out[0] ^ out[2];
   hl2[1] = hl2[1] ^ out[1] ^ out[3];

   out[1] ^= hl2[0];
   out[2] ^= hl2[1];
}

static uint32_t radix_16_truncate_asm(uint32_t a, uint32_t b) {
    uint32_t mask32 = 0x11111111;
    uint64_t a0 = (a)&mask32;
    uint64_t b0 = (b)&mask32;
    uint64_t b1 = (b >> 2) & mask32;
    uint64_t b2 = (b >> 1) & mask32;
    uint64_t b3 = (b >> 3) & mask32;
    uint32_t tmp[5] = {0};
    uint32_t out = 0;
    // 0, [out[0], tmp[1]] = a0 * b0
    _store(&out, &tmp[1], a0 * b0);
    out &= mask32;
    // 1, [tmp[2], tmp[3]] = a0 * b2
    _store(&tmp[2], &tmp[3], a0 * b2);
    // 2, [tmp[1], tmp[2]] += a0 * b1
    _store(&tmp[1], &tmp[2], _load(tmp[1], tmp[2]) + a0 * b1);
    tmp[1] &= mask32;
    uint64_t a1 = (a >> 2) & mask32;
    // 3, [tmp[1], tmp[2]] += a1 * b0
    _store(&tmp[1], &tmp[2], _load(tmp[1], tmp[2]) + a1 * b0);
    tmp[1] &= mask32;
    tmp[2] &= mask32;
    out ^= (tmp[1]<<2);
    // 4, [tmp[3], tmp[4]] += a0 * b3
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a0 * b3);
    // 5, [tmp[2], tmp[3]] += a1 * b1
    _store(&tmp[2], &tmp[3], _load(tmp[2], tmp[3]) + a1 * b1);
    tmp[2] &= mask32;
    // 6, [tmp[3], tmp[4]] += a1 * b2
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a1 * b2);
    tmp[3] &= mask32;
    uint64_t a2 = (a >> 1) & mask32;
    // 7, [tmp[2], tmp[3]] += a2 * b0
    _store(&tmp[2], &tmp[3], _load(tmp[2], tmp[3]) + a2 * b0);
    tmp[2] &= mask32;
    out ^= (tmp[2]<<1);
    // 8, [tmp[3], tmp[4]] += a2 * b1
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a2 * b1);
    uint64_t a3 = (a >> 3) & mask32;
    // 9, [tmp[3], tmp[4]] += a3 * b0
    _store(&tmp[3], &tmp[4], _load(tmp[3], tmp[4]) + a3 * b0);
    tmp[3] &= mask32;
    out ^= (tmp[3]<<3);
    return out;
}

#endif


static inline void swap32(uint32_t *res, const uint32_t *x) {
    uint32_t tmp;
	SWAPMOVE_CONST(res[0], x[0], 0x0000F0F0u, 12, tmp);
	SWAPMOVE(res[0], res[0], 0x00cc00ccu, 6, tmp);
	SWAPMOVE(res[0], res[0], 0x0a0a0a0au, 3, tmp);

}
void poly_to_rd16(uint32_t *res, const uint32_t * x, uint32_t len) {
	for(uint32_t i=0;i<len;i++){
		swap32(&res[i], &x[i]);
	}
}

static inline void inv_swap32(uint32_t * x) {
    uint32_t tmp;
	SWAPMOVE(x[0], x[0], 0x0a0a0a0au, 3, tmp);
	SWAPMOVE(x[0], x[0], 0x00cc00ccu, 6, tmp);
	SWAPMOVE(x[0], x[0], 0x0000F0F0u, 12, tmp);
}
void rd16_to_poly(uint32_t * x, uint32_t len) {
	for(uint32_t i=0;i<len;i++){
		inv_swap32(&x[i]);
	}
}

static inline void mul2_truncate (uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hl1, hl2;

    radix_16_basemul_asm(c, a[0], b[0]);
    hl1 = radix_16_truncate_asm(a[1], b[0]);
    hl2 = radix_16_truncate_asm(a[0], b[1]);

    c[1] = c[1] ^ hl1 ^ hl2;
}

static inline void mul4(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
   uint32_t hs0[2], hs1[2];
   uint32_t hl2[4];

   hs0[0] = a[0] ^ a[2];
   hs0[1] = a[1] ^ a[3];
   hs1[0] = b[0] ^ b[2];
   hs1[1] = b[1] ^ b[3];

   mul2_asm(c, a, b);
   mul2_asm(c+4, a+2, b+2);
   mul2_asm(hl2, hs0, hs1);

   hl2[0] = hl2[0] ^ c[0] ^ c[4];
   hl2[1] = hl2[1] ^ c[1] ^ c[5];
   hl2[2] = hl2[2] ^ c[2] ^ c[6];
   hl2[3] = hl2[3] ^ c[3] ^ c[7];

   c[2] ^= hl2[0];
   c[3] ^= hl2[1];
   c[4] ^= hl2[2];
   c[5] ^= hl2[3];
}

static inline void mul4_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hl1[2];
    uint32_t hl2[2];

    mul2_asm(c, a, b);
    mul2_truncate(hl1, a, b+2);
    mul2_truncate(hl2, a+2, b);

    c[2] = c[2] ^ hl1[0] ^ hl2[0];
    c[3] = c[3] ^ hl1[1] ^ hl2[1];

}

static inline void mul8(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
   uint32_t hs0[4], hs1[4];
   uint32_t hl2[8];

   hs0[0] = a[0] ^ a[4];
   hs0[1] = a[1] ^ a[5];
   hs0[2] = a[2] ^ a[6];
   hs0[3] = a[3] ^ a[7];
   hs1[0] = b[0] ^ b[4];
   hs1[1] = b[1] ^ b[5];
   hs1[2] = b[2] ^ b[6];
   hs1[3] = b[3] ^ b[7];

   mul4(c, a, b);
   mul4(c+8, a+4, b+4);
   mul4(hl2, hs0, hs1);

   hl2[0] = hl2[0] ^ c[0] ^ c[8];
   hl2[1] = hl2[1] ^ c[1] ^ c[9];
   hl2[2] = hl2[2] ^ c[2] ^ c[10];
   hl2[3] = hl2[3] ^ c[3] ^ c[11];
   hl2[4] = hl2[4] ^ c[4] ^ c[12];
   hl2[5] = hl2[5] ^ c[5] ^ c[13];
   hl2[6] = hl2[6] ^ c[6] ^ c[14];
   hl2[7] = hl2[7] ^ c[7] ^ c[15];

   c[4]  ^= hl2[0];
   c[5]  ^= hl2[1];
   c[6]  ^= hl2[2];
   c[7]  ^= hl2[3];
   c[8]  ^= hl2[4];
   c[9]  ^= hl2[5];
   c[10] ^= hl2[6];
   c[11] ^= hl2[7];
// warning #13200: No EMMS instruction before return
}

static inline void mul8_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hl1[4];
    uint32_t hl2[4];

    mul4(c, a, b);
    mul4_truncate(hl1, a, b+4);
    mul4_truncate(hl2, a+4, b);

    c[4] = c[4] ^ hl1[0] ^ hl2[0];
    c[5] = c[5] ^ hl1[1] ^ hl2[1];
    c[6] = c[6] ^ hl1[2] ^ hl2[2];
    c[7] = c[7] ^ hl1[3] ^ hl2[3];

}

// 512 x 512 -> 1024
static inline void mul16(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
   uint32_t hs0[8], hs1[8];
   uint32_t hl2[16];
    for(int i=0;i<8;i++){
        hs0[i] = a[i] ^ a[i+8];
        hs1[i] = b[i] ^ b[i+8];
    }
    mul8(c, a, b);
    mul8(c+16, a+8, b+8);
    mul8(hl2, hs0, hs1);
    for(int i=0;i<16;i++){
        hl2[i] = hl2[i] ^ c[i] ^ c[i+16];
    }
    for(int i=0;i<16;i++){
        c[i+8] ^= hl2[i];
    }
}

static inline void mul16_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hl1[8];
    uint32_t hl2[8];

    mul8(c, a, b);
    mul8_truncate(hl1, a, b+8);
    mul8_truncate(hl2, a+8, b);

    for(int i = 0; i < 8; i++)
    {
        c[i+8] = c[i+8] ^ hl1[i] ^ hl2[i];
    }

}

// 1024 x 1024 -> 2048
static inline void mul32(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
   uint32_t hs0[16], hs1[16];
    uint32_t hl2[32];
    for(int i=0;i<16;i++){
        hs0[i] = a[i] ^ a[i+16];
        hs1[i] = b[i] ^ b[i+16];
    }
    mul16(c, a, b);
    mul16(c+32, a+16, b+16);
    mul16(hl2, hs0, hs1);
    for(int i=0;i<32;i++){
        hl2[i] = hl2[i] ^ c[i] ^ c[i+32];
    }
    for(int i=0;i<32;i++){
        c[i+16] ^= hl2[i];
    }
}

static inline void mul32_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hl1[16];
    uint32_t hl2[16];

    mul16(c, a, b);
    mul16_truncate(hl1, a, b+16);
    mul16_truncate(hl2, a+16, b);

    for(int i = 0; i < 16; i++)
    {
        c[i+16] = c[i+16] ^ hl1[i] ^ hl2[i];
    }

}

// 2048 x 2048 -> 2048
static inline void mul64_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hs0[32], hs1[32];
    mul32(c, a, b);
    mul32_truncate(hs0, a+32, b);
    mul32_truncate(hs1, a, b+32);
    for(int i=0;i<32;i++){
        c[i+32] = c[i+32] ^ hs0[i] ^ hs1[i];
    }
}

// 2048 x 2048 -> 4096
static inline void mul64(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hs0[32], hs1[32];
    uint32_t hl2[64];
    for(int i=0;i<32;i++){
        hs0[i] = a[i] ^ a[i+32];
        hs1[i] = b[i] ^ b[i+32];
    }
    mul32(c, a, b);
    mul32(c+64, a+32, b+32);
    mul32(hl2, hs0, hs1);
    for(int i=0;i<64;i++){
        hl2[i] = hl2[i] ^ c[i] ^ c[i+64];
    }
    for(int i=0;i<64;i++){
        c[i+32] ^= hl2[i];
    }
}

// 4096 x 4096 -> 4096
static inline void mul128_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hs0[64], hs1[64];
    mul64(c, a, b);
    mul64_truncate(hs0, a+64, b);
    mul64_truncate(hs1, a, b+64);
    for(int i=0;i<64;i++){
        c[i+64] = c[i+64] ^ hs0[i] ^ hs1[i];
    }
}

// 4096 x 4096 -> 8192
static inline void mul128(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hs0[64], hs1[64];
    uint32_t hl2[128];
    for(int i=0;i<64;i++){
        hs0[i] = a[i] ^ a[i+64];
        hs1[i] = b[i] ^ b[i+64];
    }
    mul64(c, a, b);
    mul64(c+128, a+64, b+64);
    mul64(hl2, hs0, hs1);
    for(int i=0;i<128;i++){
        hl2[i] = hl2[i] ^ c[i] ^ c[i+128];
    }
    for(int i=0;i<128;i++){
        c[i+64] ^= hl2[i];
    }
}

// 8192 x 8192 -> 8192
static inline void mul256_truncate(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t hs0[128], hs1[128];
    mul128(c, a, b);
    mul128_truncate(hs0, a+128, b);
    mul128_truncate(hs1, a, b+128);
    for(int i=0;i<128;i++){
        c[i+128] = c[i+128] ^ hs0[i] ^ hs1[i];
    }
}


void gf2x_mul_8192(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    // This extra buffer stk[256] (and thus the memcpy in the end) is
    // necessary, because mul256_truncate and its subroutines break if c[256]
    // overlaps with a[256] or b[256].
    uint32_t stk[256];
    uint32_t sa[256], sb[256];
    
	poly_to_rd16(sa, a, 256);
	poly_to_rd16(sb, b, 256);
	mul256_truncate(stk, sa, sb);
	rd16_to_poly(stk, 256);
    memcpy(c, stk, 256 * sizeof(uint32_t));
}

void gf2x_mul_8192_no_inv_swap(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    uint32_t sa[256], sb[256];
    
	poly_to_rd16(sa, a, 256);
	poly_to_rd16(sb, b, 256);
	mul256_truncate(c, sa, sb);
}

void gf2x_mul_8192_preswap(uint32_t *c, const uint32_t *a, const uint32_t *b)
{
    // This extra buffer stk[256] (and thus the memcpy in the end) is
    // necessary, because mul256_truncate and its subroutines break if c[256]
    // overlaps with a[256] or b[256].
    uint32_t stk[256];
    mul256_truncate(stk, a, b);
    rd16_to_poly(stk, 256);
    memcpy(c, stk, 256 * sizeof(uint32_t));
}
