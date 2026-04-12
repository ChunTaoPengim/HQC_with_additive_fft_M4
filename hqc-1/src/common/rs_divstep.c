#include "fft.h"
#include "gf.h"
#include "gf256.h"
#include "parameters.h"
#include <stdint.h>
#include <string.h>
#include "rs_divstep.h"

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

static inline void swap32(uint32_t *res, const uint32_t *x) {
    uint32_t tmp;
	SWAPMOVE_CONST(res[0], x[0], 0x0000F0F0u, 12, tmp);
	SWAPMOVE(res[0], res[0], 0x00cc00ccu, 6, tmp);
	SWAPMOVE(res[0], res[0], 0x0a0a0a0au, 3, tmp);

}
static inline void inv_swap32(uint32_t * x) {
    uint32_t tmp;
	SWAPMOVE(x[0], x[0], 0x0a0a0a0au, 3, tmp);
	SWAPMOVE(x[0], x[0], 0x00cc00ccu, 6, tmp);
	SWAPMOVE(x[0], x[0], 0x0000F0F0u, 12, tmp);
}

uint32_t radix_16_gf256_4x4(uint32_t a, uint32_t b);

void compute_elp_divstep(uint8_t * omega, uint8_t *sigma,  uint8_t *syndromes)
{
    // reverse
    uint8_t g[2*(PARAM_DELTA+1)] = {0};
    memcpy(g, syndromes, 2*PARAM_DELTA);
    // 
    uint8_t v[2*(PARAM_DELTA+1)] = {0};
    uint8_t r[2*(PARAM_DELTA+1)] = {0};
    uint8_t f[2*(PARAM_DELTA+1)] = {0};
    r[0] = 1;
    f[(2*PARAM_DELTA)] = 1;
    int delta = 0;
    for (int i = 0; i < 2*PARAM_DELTA; i++) {
        // g = g*x, r = r*x
        // for(int j = 2*(PARAM_DELTA+1)-1; j > 0; j--) {
        //     g[j] = g[j-1];
        //     r[j] = r[j-1];
        // }
        memmove(g+1, g, (2*PARAM_DELTA));
        memmove(r+1, r, (2*PARAM_DELTA));
        g[0] = 0;
        r[0] = 0;
        delta -= 1;
        
        uint16_t mask_swap = ((uint16_t)(-( (g[(2*PARAM_DELTA)] != 0) && (delta < 0) ))) ;
        // conditional swap f, g
        uint8_t temp;
        for(int j = 0; j < 2*(PARAM_DELTA) +1; j++)
        {
            temp = f[j] ^ g[j];
            temp &= (uint8_t)mask_swap;
            f[j] ^= temp;
            g[j] ^= temp;
        }
        
        // conditional swap v, r
        for(int j = 0; j < 2*(PARAM_DELTA) +1; j++)
        {
            temp = v[j] ^ r[j];
            temp &= (uint8_t)mask_swap;
            r[j] ^= temp;
            v[j] ^= temp;
        }
        // delta < 0 && g[30] != 0
        delta += mask_swap & ( -delta );
        uint8_t f0 = f[2*(PARAM_DELTA)];
        uint8_t g0 = g[2*(PARAM_DELTA)];

        // g = (f2t*g + g2t*f)
        // for(int j = 0; j < 2*(PARAM_DELTA) +1; j++)
        // {
        //     g[j] = PQCLEAN_HQC128_CLEAN_gf_mul(g[j], f0) ^ PQCLEAN_HQC128_CLEAN_gf_mul(f[j], g0);
        // }
        gf256_mul((uint32_t*)g, (uint32_t)f0);
        gf256_madd((uint32_t*)g, (uint32_t*)f, (uint32_t)g0);
        
        
        // r = (f2t*r + g2t*v)
        // for(int j = 0; j < 2*(PARAM_DELTA) +1; j++)
        // {
        //     r[j] = PQCLEAN_HQC128_CLEAN_gf_mul(r[j], f0) ^ PQCLEAN_HQC128_CLEAN_gf_mul(v[j], g0);
        // }
        gf256_mul((uint32_t*)r, (uint32_t)f0);
        gf256_madd((uint32_t*)r, (uint32_t*)v, (uint32_t)g0);
        
    }
    // int degree = 0;
    // for(int i = 0; i < 32; i++)
    // {
    //     if(r[i] != 0)
    //     {
    //         degree ++;
    //     }
    // }
    // for(int i = 0; i < degree; i++)
    // {
    //     temp = r[2*PARAM_DELTA +1 - degree +i ];
    //     r[2*PARAM_DELTA +1 - degree +i ] = r[PARAM_DELTA + i];
    //     r[PARAM_DELTA + i] = temp; 
    // }
    
    memcpy(sigma, r+PARAM_DELTA, PARAM_DELTA + 1);
    memcpy(omega, g+PARAM_DELTA, (PARAM_DELTA+1));
}

static void poly_eval(uint8_t * sigma, uint8_t * omega, uint8_t * beta_j_inverse, uint8_t * tmpp1_array_bytes, uint8_t * tmpp2_array_bytes)
{
    uint32_t tmpp1_array[(PARAM_DELTA + 3)/4] = {0};
    uint32_t tmpp2_array[(PARAM_DELTA + 3)/4] = {0};

    uint32_t * beta_j_inverse_u32 = (uint32_t *)beta_j_inverse;
    
    gf256_madd_16(tmpp1_array, beta_j_inverse_u32, omega[(PARAM_DELTA-1)]);
    gf256_madd_16(tmpp2_array, beta_j_inverse_u32, sigma[(PARAM_DELTA-1)]);
    tmpp1_array[0] ^= (uint32_t)omega[(PARAM_DELTA-2)] * 0x01010101u;
    tmpp2_array[0] ^= (uint32_t)sigma[(PARAM_DELTA-2)] * 0x01010101u;


    // pre-compute (uint32_t)omega[j] * 0x01010101u
    // pre-compute (uint32_t)sigma[j] * 0x01010101u
    uint32_t omegax01010101[PARAM_DELTA-2] = {0};
    uint32_t sigmax01010101[PARAM_DELTA-2] = {0};
    uint32_t omegax01010101_rdx[PARAM_DELTA-2] = {0};
    uint32_t sigmax01010101_rdx[PARAM_DELTA-2] = {0};
    for(int j = PARAM_DELTA-3; j >= 0; j--){
        omegax01010101[j] = (uint32_t)omega[j] * 0x01010101u;
        swap32(&omegax01010101_rdx[j], &omegax01010101[j]);
    }
    for(int i = 0; i <= PARAM_DELTA-3; i += 2){
        sigmax01010101[i] = (uint32_t)sigma[i] * 0x01010101u;
        swap32(&sigmax01010101_rdx[i], &sigmax01010101[i]);
    }
    

    uint32_t beta_j_inv_rdx;
    for(int i = 0; i < 4; i++)
    {
        swap32(&tmpp1_array[i], &tmpp1_array[i]);
        swap32(&tmpp2_array[i], &tmpp2_array[i]);
        swap32(&beta_j_inv_rdx, &beta_j_inverse_u32[i]);
        for(int j = PARAM_DELTA-3; j >= 1; j-=2) // PARAM_DELTA = 15 in HQC-1
        {
            tmpp1_array[i] = radix_16_gf256_4x4(tmpp1_array[i], beta_j_inv_rdx);
            tmpp2_array[i] = radix_16_gf256_4x4(tmpp2_array[i], beta_j_inv_rdx);
            tmpp1_array[i] ^= omegax01010101_rdx[j];
            tmpp2_array[i] ^= sigmax01010101_rdx[j];
            tmpp1_array[i] = radix_16_gf256_4x4(tmpp1_array[i], beta_j_inv_rdx);
            tmpp2_array[i] = radix_16_gf256_4x4(tmpp2_array[i], beta_j_inv_rdx);
            tmpp1_array[i] ^= omegax01010101_rdx[j-1];
            // skip below because they're all zeroes
            //tmpp2_array[i] ^= sigmax01010101_rdx[j-1];
        }
        tmpp1_array[i] = radix_16_gf256_4x4(tmpp1_array[i], beta_j_inv_rdx);
        tmpp2_array[i] = radix_16_gf256_4x4(tmpp2_array[i], beta_j_inv_rdx);
        tmpp1_array[i] ^= omegax01010101_rdx[0];
        tmpp2_array[i] ^= sigmax01010101_rdx[0];

        inv_swap32(&tmpp1_array[i]);
        inv_swap32(&tmpp2_array[i]);
    }
    memcpy(tmpp1_array_bytes, tmpp1_array, PARAM_DELTA);
    memcpy(tmpp2_array_bytes, tmpp2_array, PARAM_DELTA);
}
void compute_error_values_new(uint8_t *error_values, const uint8_t *omega, const uint8_t *sigma, const uint8_t *error) {
    // uint8_t beta_j[PARAM_DELTA+1] = {0};
    uint8_t beta_j_inverse[PARAM_DELTA+1] = {0};
    uint8_t e_j[PARAM_DELTA + 1] = {0};

    uint16_t delta_counter;
    uint16_t found;
    uint16_t mask1;
    uint16_t mask2;
    // uint32_t tmpp1_array[(PARAM_DELTA + 3)/4] = {0};
    // uint32_t tmpp2_array[(PARAM_DELTA + 3)/4] = {0};
    uint8_t tmpp1_array[PARAM_DELTA] = {0};
    uint8_t tmpp2_array[PARAM_DELTA] = {0};

    // Compute the beta_{j_i} page 31 of the documentation
    delta_counter = 0;
    for (size_t i = 0; i < PARAM_N1; i++) {
        found = 0;
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0
        for (size_t j = 0; j < PARAM_DELTA; j++) {
            mask2 = ~((uint16_t) (-((int32_t) j ^ delta_counter) >> 31)); // j == delta_counter
            beta_j_inverse[j] += mask1 & mask2 & gf_exp_inverse[i];
            found += mask1 & mask2 & 1;
        }
        delta_counter += found;
    }
    

    // compute the derivative of sigma
    uint8_t sigma_deriv[PARAM_DELTA + 1] = {0};
    for(int i = 1; i <= PARAM_DELTA; i += 2) {
        sigma_deriv[i - 1] = sigma[i];
    }
    poly_eval((uint8_t *)sigma_deriv, (uint8_t *)omega, beta_j_inverse, tmpp1_array, tmpp2_array);
    for(int i = 0; i < PARAM_DELTA; i++)
    {
        e_j[i] = PQCLEAN_HQC128_CLEAN_gf_mul(tmpp1_array[i], PQCLEAN_HQC128_CLEAN_gf_inverse(tmpp2_array[i]));
    }
    // 

    delta_counter = 0;
    for (size_t i = 0; i < PARAM_N1; ++i) {
        found = 0;
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0
        for (size_t j = 0; j < PARAM_DELTA; j++) {
            mask2 = ~((uint16_t) (-((int32_t) j ^ delta_counter) >> 31)); // j == delta_counter
            error_values[i] += mask1 & mask2 & e_j[j];
            found += mask1 & mask2 & 1;
        }
        delta_counter += found;
    }
}