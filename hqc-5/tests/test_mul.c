#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "hal.h"
#include "api.h"
#include "parameters.h"
#include "gfmul_fft.h"
#include "gf2x.h"
#include "reduce.h"


static void schoolbook_mul(uint32_t *r, const uint32_t *a, const uint32_t *b, size_t n) {
    memset(r, 0, 2 * n * sizeof(uint32_t));
    for (size_t i = 0; i < n; i++) {
        uint32_t ai = a[i];
        for (int bit = 0; bit < 32; bit++) {
            uint32_t mask = -((ai >> bit) & 1U);
            size_t base = i;
            int sh = bit;
            int inv = 32 - sh;
            if (sh == 0) {
                for (size_t j = 0; j < n; j++) {
                    r[base + j] ^= b[j] & mask;
                }
            } else {
                for (size_t j = 0; j < n; j++) {
                    r[base + j] ^= (b[j] << sh) & mask;
                    r[base + j + 1] ^= (b[j] >> inv) & mask;
                }
            }
        }
    }
}

/*
testing between schoolbook multiplication and FFT multiplication
*/
void test1()
{
    uint64_t y[VEC_N_SIZE_64] = {0};
    uint64_t h[VEC_N_SIZE_64] = {0};
    uint64_t s[VEC_N_SIZE_64] = {0};
    for(int i = 0; i < VEC_N_SIZE_64 - 1; i++) {
        y[i] = rand();
        h[i] = rand();
    }
    PQCLEAN_HQC256_CLEAN_vect_mul(s, y, h);
    uint64_t s_2[VEC_N_SIZE_64*2] = {0};
    // schoolbook multiplication
    schoolbook_mul((uint32_t *)s_2, (uint32_t *)y, (uint32_t *)h, VEC_N_SIZE_64 * 2);
    // reduction
    uint64_t s_reduced[VEC_N_SIZE_64] = {0};
    reduce((uint32_t*)s_reduced, (uint32_t *)s_2);
    // compare
    for(int i = 0; i < VEC_N_SIZE_64; i++)
    {
        if(s[i] != s_reduced[i]) {
            char outstr[128];
            sprintf(outstr, "Error at index %d: s=%llu, s_reduced=%llu\n", i, s[i], ((uint64_t *)s_reduced)[i]);
            hal_send_str(outstr);
        }
    }
    hal_send_str("\n*** ALL GOOD ***\n");

}

int main()
{
    hal_setup(CLOCK_BENCHMARK);
    test1();
    return 0;
}
