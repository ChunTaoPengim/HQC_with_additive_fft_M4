#ifndef _GF256_RDX16_H_
#define _GF256_RDX16_H_

#include <stdint.h>

#include "parameters.h" // for defining _M4_ASM_

#define _0x1D_RDX16_  0x00011101
#define _NORM_RDX16_  0x11111111
#define _U64NORM_RDX16_  0x1111111111111111


#if defined(_M4_ASM_)

static inline
uint64_t rdx16_mul( uint32_t a , uint32_t b ) {
    uint32_t lo, hi;
    __asm__ volatile (
    "umull %0, %1, %2, %3"
    : "=r" (lo), "=r" (hi)
    : "r" (a), "r" (b)
    );
    return ((uint64_t)hi << 32) | lo;
}

static inline
uint64_t rdx16_muladd( uint64_t c , uint32_t a , uint32_t b ) {
    uint32_t lo = (uint32_t)(c & 0xFFFFFFFF);
    uint32_t hi = (uint32_t)(c >> 32);
    __asm__ volatile (
    "umlal %0, %1, %2, %3"
    : "+r" (lo), "+r" (hi)
    : "r" (a), "r" (b)
    );
    return ((uint64_t)hi << 32) | lo;
}

static inline
uint32_t mla( uint32_t a , uint32_t b , uint32_t c )
{
    uint32_t result;
    __asm__ volatile (
        "mla %0, %1, %2, %3"
        : "=r" (result)
        : "r" (a),
          "r" (b),
          "r" (c)
    );
    return result;
}

#else

static inline
uint64_t rdx16_mul( uint32_t a , uint32_t b ) { return ((uint64_t)(a)) * ((uint64_t)(b)); }

static inline
uint64_t rdx16_muladd( uint64_t c , uint32_t a , uint32_t b ) { return c + (((uint64_t)(a)) * ((uint64_t)(b))); }

static inline
uint32_t mla( uint32_t a , uint32_t b , uint32_t c ) { return (a*b)+c; }

#endif

// 0x11d GF(256) reduction for rdx16 representation
static inline
uint32_t rdx16_gf256_reduce( uint32_t x_hi8bit )  // presume: x_hi8bit contain only 7-bit nomialize rdx16 input
{
    uint64_t r1 = rdx16_mul( x_hi8bit , _0x1D_RDX16_ );
    uint32_t r1h = (r1>>32);  // high 32-bit register
    uint32_t r1l = (r1&0xFFFFFFFF); // low 32-bit register
    uint32_t r2l = mla( r1h, _0x1D_RDX16_ , r1l );  // multiplication without normalization
    return r2l;   // no nomailize in the result.
}

static inline
uint64_t rdx16_gf256_mul_without_reduce( uint32_t a , uint32_t b ) { return rdx16_mul( a , b )& _U64NORM_RDX16_; }




///////////////////////  data conversion ///////////////////////



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

#endif


static inline
uint32_t rdx16_from_bitseq( uint32_t a )
{
    uint32_t res;
    uint32_t tmp;
	SWAPMOVE_CONST(res, a, 0x0000F0F0u, 12, tmp);
	SWAPMOVE(res, res, 0x00cc00ccu, 6, tmp);
	SWAPMOVE(res, res, 0x0a0a0a0au, 3, tmp);
    return res;
}

static inline
uint32_t rdx16_to_bitseq( uint32_t a )
{
    uint32_t tmp;
	SWAPMOVE(a, a, 0x0a0a0a0au, 3, tmp);
	SWAPMOVE(a, a, 0x00cc00ccu, 6, tmp);
	SWAPMOVE(a, a, 0x0000F0F0u, 12, tmp);
    return a;
}







#endif /* _GF256_RDX16_H_ */

