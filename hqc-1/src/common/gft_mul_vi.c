// Implemented by Myeonghoon Lee and Jihoon Jang
// public domain

#include "parameters.h"

#include "gft_mul_vi.h"
#include "gen/gft_mul_s_inv_i_beta_matrix_defs.h"

#if !defined(_M4_ASM_)

static void mul_matrix( sto_t * out , const sto_t * in , const uint8_t matrix[32][32] ) {
    for (int i = 0; i < 32; i++) {
        out[i] = 0;
        for (int j = 0; j < 32; j++) {
            out[i] ^= matrix[i][j] * in[j];
        }
    }
}

#define MUL_S_INV_I_BETA(IDX)                                           \
    void gft_mul_s_inv_##IDX##_beta( sto_t * out , const sto_t * in ) { \
        mul_matrix(out, in, gft_mul_s_inv_##IDX##_beta_matrix);         \
    }

MUL_S_INV_I_BETA(1)
MUL_S_INV_I_BETA(2)
MUL_S_INV_I_BETA(3)
MUL_S_INV_I_BETA(4)
MUL_S_INV_I_BETA(5)
MUL_S_INV_I_BETA(6)
MUL_S_INV_I_BETA(7)
MUL_S_INV_I_BETA(8)
MUL_S_INV_I_BETA(9)
MUL_S_INV_I_BETA(10)
MUL_S_INV_I_BETA(11)

#endif  //!defined(_M4_ASM_)
