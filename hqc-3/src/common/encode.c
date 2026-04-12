// Implemented by Ming-Shing Chen, Tung Chou and Markus Krausz.
// Modified by Till Eifert
// public domain

#include "encode.h"

#include "parameters.h"
#include "gen/encode_matrix_defs.h"

#if !defined(_M4_ASM_)
void encode_to_gft_full_length(uint32_t * out , const uint32_t * in ) {
    for (int i = 0; i < 32; i++) {
        out[i] = 0;
        for (int j = 0; j < 32; j++) {
            out[i] ^= encode_matrix[i][j] * in[j];
        }
    }
}

void decode_from_gft_full_length( uint32_t * out , const uint32_t * in ) {
    for (int i = 0; i < 32; i++) {
        out[i] = 0;
        for (int j = 0; j < 32; j++) {
            out[i] ^= decode_matrix[i][j] * in[j];
        }
    }
}
#endif
