#include "fft.h"
#include "gf.h"
#include "parameters.h"
#include "reed_solomon.h"
#include "gf256_rdx16.h"
#include "gf256.h"
#include <stdint.h>
#include <string.h>
/**
 * @file reed_solomon.c
 * @brief Constant time implementation of Reed-Solomon codes
 */

#if 1

#include "rs_encoding.h"
#include "rs_divstep.h"

/**
 * @brief Encodes a message message of PARAM_K bits to a Reed-Solomon codeword codeword of PARAM_N1 bytes
 *
 * Following @cite lin1983error (Chapter 4 - Cyclic Codes),
 * We perform a systematic encoding using a linear (PARAM_N1 - PARAM_K)-stage shift register
 * with feedback connections based on the generator polynomial PARAM_RS_POLY of the Reed-Solomon code.
 *
 * @param[out] cdw Array of size VEC_N1_SIZE_64 receiving the encoded message
 * @param[in] msg Array of size VEC_K_SIZE_64 storing the message
 */
void PQCLEAN_HQC128_CLEAN_reed_solomon_encode(uint8_t *cdw, const uint8_t *msg) {
#if (46==PARAM_N1)&&(16==PARAM_K)
    rs_encode_46_16(cdw,msg);
#elif (56==PARAM_N1)&&(24==PARAM_K)
    rs_encode_56_24(cdw,msg);
#elif (90==PARAM_N1)&&(32==PARAM_K)
    rs_encode_90_32(cdw,msg);
#else
error here: no supported encoding parameter.
#endif
}

#else

/**
 * @brief Encodes a message message of PARAM_K bits to a Reed-Solomon codeword codeword of PARAM_N1 bytes
 *
 * Following @cite lin1983error (Chapter 4 - Cyclic Codes),
 * We perform a systematic encoding using a linear (PARAM_N1 - PARAM_K)-stage shift register
 * with feedback connections based on the generator polynomial PARAM_RS_POLY of the Reed-Solomon code.
 *
 * @param[out] cdw Array of size VEC_N1_SIZE_64 receiving the encoded message
 * @param[in] msg Array of size VEC_K_SIZE_64 storing the message
 */
void PQCLEAN_HQC128_CLEAN_reed_solomon_encode(uint8_t *cdw, const uint8_t *msg) {
    uint8_t gate_value = 0;

    uint16_t tmp[PARAM_G] = {0};
    uint16_t PARAM_RS_POLY [] = {RS_POLY_COEFS};

    memset(cdw, 0, PARAM_N1);

    for (size_t i = 0; i < PARAM_K; ++i) {
        gate_value = msg[PARAM_K - 1 - i] ^ cdw[PARAM_N1 - PARAM_K - 1];

        for (size_t j = 0; j < PARAM_G; ++j) {
            tmp[j] = PQCLEAN_HQC128_CLEAN_gf_mul(gate_value, PARAM_RS_POLY[j]);
        }

        for (size_t k = PARAM_N1 - PARAM_K - 1; k; --k) {
            cdw[k] = (uint8_t)(cdw[k - 1] ^ tmp[k]);
        }

        cdw[0] = (uint8_t)tmp[0];
    }
    memcpy(cdw + PARAM_N1 - PARAM_K, msg, PARAM_K);
}

#endif


#if 1

#include "rs_syndrome.h"

/**
 * @brief Computes 2 * PARAM_DELTA syndromes
 *
 * @param[out] syndromes Array of size 2 * PARAM_DELTA receiving the computed syndromes
 * @param[in] cdw Array of size PARAM_N1 storing the received vector
 */
static void compute_syndromes(uint8_t *syndromes, uint8_t *cdw) {
#if (46==PARAM_N1)&&(16==PARAM_K)
    rs_syndrome_n46_r30(syndromes,cdw);
#elif (56==PARAM_N1)&&(24==PARAM_K)
    rs_syndrome_n56_r32(syndromes,cdw);
#elif (90==PARAM_N1)&&(32==PARAM_K)
    rs_syndrome_n90_r58(syndromes,cdw);
#else
error here: no supported encoding parameter.
#endif

}

#elif 0

static const uint8_t alpha_ij_pow_trans[45][32] ={
    {2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 0, 0},
    {4, 16, 64, 29, 116, 205, 19, 76, 45, 180, 234, 143, 6, 24, 96, 157, 78, 37, 148, 106, 181, 238, 159, 70, 5, 20, 80, 93, 105, 185, 0, 0},
    {8, 64, 58, 205, 38, 45, 117, 143, 12, 96, 39, 37, 53, 181, 193, 70, 10, 80, 186, 185, 161, 97, 47, 101, 15, 120, 231, 107, 127, 223, 0, 0},
    {16, 29, 205, 76, 180, 143, 24, 157, 37, 106, 238, 70, 20, 93, 185, 95, 153, 101, 30, 253, 107, 254, 91, 217, 17, 13, 208, 129, 248, 59, 0, 0},
    {32, 116, 38, 180, 3, 96, 156, 106, 193, 5, 160, 185, 190, 94, 15, 253, 214, 223, 226, 17, 26, 103, 124, 59, 51, 46, 169, 132, 77, 85, 0, 0},
    {64, 205, 45, 143, 96, 37, 181, 70, 80, 185, 97, 101, 120, 107, 223, 217, 68, 208, 62, 59, 102, 184, 33, 168, 85, 228, 191, 252, 241, 150, 0, 0},
    {128, 19, 117, 24, 156, 181, 140, 93, 161, 94, 60, 107, 163, 67, 26, 129, 147, 102, 109, 132, 41, 57, 209, 252, 255, 98, 87, 200, 224, 89, 0, 0},
    {29, 76, 143, 157, 106, 70, 93, 95, 101, 253, 254, 217, 13, 129, 59, 133, 79, 168, 73, 230, 252, 227, 149, 130, 28, 81, 195, 18, 247, 44, 0, 0},
    {58, 45, 12, 37, 193, 80, 161, 101, 231, 223, 134, 208, 237, 102, 169, 168, 146, 191, 179, 150, 87, 7, 166, 195, 36, 251, 125, 173, 64, 38, 0, 0},
    {116, 180, 96, 106, 5, 185, 94, 253, 223, 17, 103, 59, 46, 132, 85, 230, 215, 150, 174, 28, 89, 172, 244, 44, 108, 32, 38, 3, 156, 193, 0, 0},
    {232, 234, 39, 238, 160, 97, 60, 254, 134, 103, 118, 184, 84, 57, 145, 227, 220, 7, 162, 172, 245, 176, 71, 58, 180, 192, 181, 40, 95, 15, 0, 0},
    {205, 143, 37, 70, 185, 101, 107, 217, 208, 59, 184, 168, 228, 252, 150, 130, 221, 195, 61, 44, 173, 58, 117, 39, 193, 186, 47, 231, 182, 26, 0, 0},
    {135, 6, 53, 20, 190, 120, 163, 13, 237, 46, 84, 228, 229, 98, 100, 81, 69, 251, 131, 32, 45, 192, 238, 186, 94, 187, 217, 189, 236, 169, 0, 0},
    {19, 24, 181, 93, 94, 107, 67, 129, 102, 132, 57, 252, 98, 200, 89, 18, 11, 173, 232, 3, 53, 40, 194, 231, 226, 189, 197, 158, 170, 145, 0, 0},
    {38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 0, 0},
    {76, 157, 70, 95, 253, 217, 129, 133, 168, 230, 227, 130, 81, 18, 44, 2, 152, 39, 140, 190, 231, 175, 31, 23, 77, 209, 219, 25, 162, 36, 0, 0},
    {152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 0, 0},
    {45, 37, 80, 101, 223, 208, 102, 168, 191, 150, 7, 195, 251, 173, 38, 39, 10, 47, 127, 26, 197, 21, 115, 219, 100, 242, 245, 54, 205, 96, 0, 0},
    {90, 148, 186, 30, 226, 62, 109, 73, 179, 174, 162, 61, 131, 232, 96, 140, 153, 127, 52, 51, 168, 99, 98, 56, 172, 22, 8, 234, 212, 185, 0, 0},
    {180, 106, 185, 253, 17, 59, 132, 230, 150, 28, 172, 44, 32, 3, 193, 190, 214, 26, 51, 77, 145, 55, 167, 36, 233, 116, 96, 5, 94, 223, 0, 0},
    {117, 181, 161, 107, 26, 102, 41, 252, 87, 89, 245, 173, 45, 53, 185, 231, 68, 197, 168, 145, 110, 166, 61, 54, 38, 37, 186, 120, 134, 59, 0, 0},
    {234, 238, 97, 254, 103, 184, 57, 227, 7, 172, 176, 58, 192, 40, 15, 175, 147, 21, 99, 55, 166, 122, 216, 45, 106, 222, 107, 52, 133, 85, 0, 0},
    {201, 159, 47, 91, 124, 33, 209, 149, 166, 244, 71, 117, 238, 194, 223, 31, 79, 115, 98, 167, 61, 216, 90, 181, 190, 254, 206, 218, 213, 150, 0, 0},
    {143, 70, 101, 217, 59, 168, 252, 130, 195, 44, 58, 39, 186, 231, 26, 23, 146, 219, 56, 36, 54, 45, 181, 97, 223, 62, 33, 191, 110, 89, 0, 0},
    {3, 5, 15, 17, 51, 85, 255, 28, 36, 108, 180, 193, 94, 226, 59, 77, 215, 100, 172, 233, 38, 106, 190, 223, 124, 132, 145, 174, 239, 44, 0, 0},
    {6, 20, 120, 13, 46, 228, 98, 81, 251, 32, 192, 186, 187, 189, 169, 209, 220, 242, 22, 116, 37, 222, 254, 62, 132, 63, 130, 43, 250, 38, 0, 0},
    {12, 80, 231, 208, 169, 191, 87, 195, 125, 38, 181, 47, 217, 197, 85, 219, 221, 245, 8, 96, 186, 107, 206, 33, 145, 130, 86, 207, 45, 193, 0, 0},
    {24, 93, 107, 129, 132, 252, 200, 18, 173, 3, 40, 231, 189, 158, 145, 25, 69, 54, 234, 5, 120, 52, 218, 191, 174, 43, 207, 90, 35, 15, 0, 0},
    {48, 105, 127, 248, 77, 241, 224, 247, 64, 156, 95, 182, 236, 170, 150, 162, 11, 205, 212, 94, 134, 133, 213, 110, 239, 250, 45, 35, 30, 26, 0, 0},
    {96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 0, 0},
    {192, 222, 182, 151, 114, 110, 155, 27, 143, 160, 177, 237, 82, 75, 89, 88, 152, 70, 240, 103, 21, 123, 224, 251, 116, 212, 101, 136, 218, 145, 0, 0},
    {157, 95, 217, 133, 230, 130, 18, 2, 39, 190, 175, 23, 209, 25, 36, 4, 78, 97, 67, 46, 191, 50, 72, 8, 156, 194, 134, 92, 99, 100, 0, 0},
    {39, 97, 134, 184, 145, 7, 245, 58, 181, 15, 208, 21, 241, 166, 44, 45, 10, 107, 237, 85, 196, 195, 54, 12, 185, 182, 102, 115, 130, 36, 0, 0},
    {78, 153, 68, 79, 215, 221, 11, 152, 10, 214, 147, 146, 220, 69, 1, 78, 153, 68, 79, 215, 221, 11, 152, 10, 214, 147, 146, 220, 69, 1, 0, 0},
    {156, 94, 26, 132, 255, 89, 233, 3, 185, 226, 46, 145, 28, 235, 38, 5, 214, 59, 114, 174, 36, 32, 106, 15, 103, 77, 150, 239, 108, 96, 0, 0},
    {37, 101, 208, 168, 150, 195, 173, 39, 47, 26, 21, 219, 242, 54, 96, 97, 68, 33, 241, 89, 207, 12, 161, 134, 169, 179, 166, 125, 143, 185, 0, 0},
    {74, 137, 206, 82, 55, 138, 16, 212, 120, 124, 73, 87, 72, 29, 193, 211, 147, 228, 25, 244, 205, 140, 177, 197, 230, 141, 251, 76, 40, 223, 0, 0},
    {148, 30, 62, 73, 174, 61, 232, 140, 127, 51, 99, 56, 22, 234, 185, 67, 79, 241, 121, 108, 39, 188, 189, 41, 55, 9, 64, 238, 211, 59, 0, 0},
    {53, 120, 237, 228, 100, 251, 45, 186, 217, 169, 241, 242, 173, 37, 15, 62, 146, 130, 245, 38, 80, 182, 184, 179, 89, 54, 39, 101, 206, 85, 0, 0},
    {106, 253, 59, 230, 28, 44, 3, 190, 26, 77, 55, 36, 116, 5, 223, 46, 215, 89, 108, 156, 15, 124, 114, 100, 235, 180, 185, 17, 132, 150, 0, 0},
    {212, 211, 197, 198, 167, 207, 157, 202, 62, 114, 200, 139, 201, 95, 26, 154, 220, 61, 19, 160, 217, 158, 171, 86, 32, 159, 127, 133, 229, 89, 0, 0},
    {181, 107, 102, 252, 89, 173, 53, 231, 197, 145, 166, 54, 37, 120, 59, 191, 221, 207, 39, 15, 237, 115, 56, 125, 96, 101, 62, 228, 7, 44, 0, 0},
    {119, 177, 23, 123, 239, 8, 159, 225, 184, 255, 43, 64, 140, 91, 169, 171, 69, 58, 20, 226, 33, 49, 18, 205, 160, 67, 21, 149, 144, 38, 0, 0},
    {238, 254, 184, 227, 172, 58, 40, 175, 21, 55, 122, 45, 222, 52, 85, 50, 11, 12, 188, 124, 115, 224, 131, 37, 253, 151, 252, 121, 2, 193, 0, 0},
    {193, 223, 169, 150, 36, 38, 185, 26, 85, 100, 44, 96, 15, 59, 145, 89, 1, 193, 223, 169, 150, 36, 38, 185, 26, 85, 100, 44, 96, 15, 0, 0}
};

/**
 * @brief Computes 2 * PARAM_DELTA syndromes
 *
 * @param[out] syndromes Array of size 2 * PARAM_DELTA receiving the computed syndromes
 * @param[in] cdw Array of size PARAM_N1 storing the received vector
 */
#include "gf256.h"
static void compute_syndromes(uint8_t *syndromes, uint8_t *cdw) {
    memset(syndromes, cdw[0], 2 * PARAM_DELTA);
    for (size_t j = 1; j < PARAM_N1; ++j) {
        gf256_madd((uint32_t*)syndromes, (uint32_t*)alpha_ij_pow_trans[j-1], cdw[j]);
    }
}

#else

static const uint16_t alpha_ij_pow [30][45] = {{2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193}, {4, 16, 64, 29, 116, 205, 19, 76, 45, 180, 234, 143, 6, 24, 96, 157, 78, 37, 148, 106, 181, 238, 159, 70, 5, 20, 80, 93, 105, 185, 222, 95, 97, 153, 94, 101, 137, 30, 120, 253, 211, 107, 177, 254, 223}, {8, 64, 58, 205, 38, 45, 117, 143, 12, 96, 39, 37, 53, 181, 193, 70, 10, 80, 186, 185, 161, 97, 47, 101, 15, 120, 231, 107, 127, 223, 182, 217, 134, 68, 26, 208, 206, 62, 237, 59, 197, 102, 23, 184, 169}, {16, 29, 205, 76, 180, 143, 24, 157, 37, 106, 238, 70, 20, 93, 185, 95, 153, 101, 30, 253, 107, 254, 91, 217, 17, 13, 208, 129, 248, 59, 151, 133, 184, 79, 132, 168, 82, 73, 228, 230, 198, 252, 123, 227, 150}, {32, 116, 38, 180, 3, 96, 156, 106, 193, 5, 160, 185, 190, 94, 15, 253, 214, 223, 226, 17, 26, 103, 124, 59, 51, 46, 169, 132, 77, 85, 114, 230, 145, 215, 255, 150, 55, 174, 100, 28, 167, 89, 239, 172, 36}, {64, 205, 45, 143, 96, 37, 181, 70, 80, 185, 97, 101, 120, 107, 223, 217, 68, 208, 62, 59, 102, 184, 33, 168, 85, 228, 191, 252, 241, 150, 110, 130, 7, 221, 89, 195, 138, 61, 251, 44, 207, 173, 8, 58, 38}, {128, 19, 117, 24, 156, 181, 140, 93, 161, 94, 60, 107, 163, 67, 26, 129, 147, 102, 109, 132, 41, 57, 209, 252, 255, 98, 87, 200, 224, 89, 155, 18, 245, 11, 233, 173, 16, 232, 45, 3, 157, 53, 159, 40, 185}, {29, 76, 143, 157, 106, 70, 93, 95, 101, 253, 254, 217, 13, 129, 59, 133, 79, 168, 73, 230, 252, 227, 149, 130, 28, 81, 195, 18, 247, 44, 27, 2, 58, 152, 3, 39, 212, 140, 186, 190, 202, 231, 225, 175, 26}, {58, 45, 12, 37, 193, 80, 161, 101, 231, 223, 134, 208, 237, 102, 169, 168, 146, 191, 179, 150, 87, 7, 166, 195, 36, 251, 125, 173, 64, 38, 143, 39, 181, 10, 185, 47, 120, 127, 217, 26, 62, 197, 184, 21, 85}, {116, 180, 96, 106, 5, 185, 94, 253, 223, 17, 103, 59, 46, 132, 85, 230, 215, 150, 174, 28, 89, 172, 244, 44, 108, 32, 38, 3, 156, 193, 160, 190, 15, 214, 226, 26, 124, 51, 169, 77, 114, 145, 255, 55, 100}, {232, 234, 39, 238, 160, 97, 60, 254, 134, 103, 118, 184, 84, 57, 145, 227, 220, 7, 162, 172, 245, 176, 71, 58, 180, 192, 181, 40, 95, 15, 177, 175, 208, 147, 46, 21, 73, 99, 241, 55, 200, 166, 43, 122, 44}, {205, 143, 37, 70, 185, 101, 107, 217, 208, 59, 184, 168, 228, 252, 150, 130, 221, 195, 61, 44, 173, 58, 117, 39, 193, 186, 47, 231, 182, 26, 237, 23, 21, 146, 145, 219, 87, 56, 242, 36, 139, 54, 64, 45, 96}, {135, 6, 53, 20, 190, 120, 163, 13, 237, 46, 84, 228, 229, 98, 100, 81, 69, 251, 131, 32, 45, 192, 238, 186, 94, 187, 217, 189, 236, 169, 82, 209, 241, 220, 28, 242, 72, 22, 173, 116, 201, 37, 140, 222, 15}, {19, 24, 181, 93, 94, 107, 67, 129, 102, 132, 57, 252, 98, 200, 89, 18, 11, 173, 232, 3, 53, 40, 194, 231, 226, 189, 197, 158, 170, 145, 75, 25, 166, 69, 235, 54, 29, 234, 37, 5, 95, 120, 91, 52, 59}, {38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145}, {76, 157, 70, 95, 253, 217, 129, 133, 168, 230, 227, 130, 81, 18, 44, 2, 152, 39, 140, 190, 231, 175, 31, 23, 77, 209, 219, 25, 162, 36, 88, 4, 45, 78, 5, 97, 211, 67, 62, 46, 154, 191, 171, 50, 89}, {152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1}, {45, 37, 80, 101, 223, 208, 102, 168, 191, 150, 7, 195, 251, 173, 38, 39, 10, 47, 127, 26, 197, 21, 115, 219, 100, 242, 245, 54, 205, 96, 70, 97, 107, 68, 59, 33, 228, 241, 130, 89, 61, 207, 58, 12, 193}, {90, 148, 186, 30, 226, 62, 109, 73, 179, 174, 162, 61, 131, 232, 96, 140, 153, 127, 52, 51, 168, 99, 98, 56, 172, 22, 8, 234, 212, 185, 240, 67, 237, 79, 114, 241, 25, 121, 245, 108, 19, 39, 20, 188, 223}, {180, 106, 185, 253, 17, 59, 132, 230, 150, 28, 172, 44, 32, 3, 193, 190, 214, 26, 51, 77, 145, 55, 167, 36, 233, 116, 96, 5, 94, 223, 103, 46, 85, 215, 174, 89, 244, 108, 38, 156, 160, 15, 226, 124, 169}, {117, 181, 161, 107, 26, 102, 41, 252, 87, 89, 245, 173, 45, 53, 185, 231, 68, 197, 168, 145, 110, 166, 61, 54, 38, 37, 186, 120, 134, 59, 21, 191, 196, 221, 36, 207, 205, 39, 80, 15, 217, 237, 33, 115, 150}, {234, 238, 97, 254, 103, 184, 57, 227, 7, 172, 176, 58, 192, 40, 15, 175, 147, 21, 99, 55, 166, 122, 216, 45, 106, 222, 107, 52, 133, 85, 123, 50, 195, 11, 32, 12, 140, 188, 182, 124, 158, 115, 49, 224, 36}, {201, 159, 47, 91, 124, 33, 209, 149, 166, 244, 71, 117, 238, 194, 223, 31, 79, 115, 98, 167, 61, 216, 90, 181, 190, 254, 206, 218, 213, 150, 224, 72, 54, 152, 106, 161, 177, 189, 184, 114, 171, 56, 18, 131, 38}, {143, 70, 101, 217, 59, 168, 252, 130, 195, 44, 58, 39, 186, 231, 26, 23, 146, 219, 56, 36, 54, 45, 181, 97, 223, 62, 33, 191, 110, 89, 251, 8, 12, 10, 15, 134, 197, 41, 179, 100, 86, 125, 205, 37, 185}, {3, 5, 15, 17, 51, 85, 255, 28, 36, 108, 180, 193, 94, 226, 59, 77, 215, 100, 172, 233, 38, 106, 190, 223, 124, 132, 145, 174, 239, 44, 116, 156, 185, 214, 103, 169, 230, 55, 89, 235, 32, 96, 160, 253, 26}, {6, 20, 120, 13, 46, 228, 98, 81, 251, 32, 192, 186, 187, 189, 169, 209, 220, 242, 22, 116, 37, 222, 254, 62, 132, 63, 130, 43, 250, 38, 212, 194, 182, 147, 77, 179, 141, 9, 54, 180, 159, 101, 67, 151, 85}, {12, 80, 231, 208, 169, 191, 87, 195, 125, 38, 181, 47, 217, 197, 85, 219, 221, 245, 8, 96, 186, 107, 206, 33, 145, 130, 86, 207, 45, 193, 101, 134, 102, 146, 150, 166, 251, 64, 39, 185, 127, 62, 21, 252, 100}, {24, 93, 107, 129, 132, 252, 200, 18, 173, 3, 40, 231, 189, 158, 145, 25, 69, 54, 234, 5, 120, 52, 218, 191, 174, 43, 207, 90, 35, 15, 136, 92, 115, 220, 239, 125, 76, 238, 101, 17, 133, 228, 149, 121, 44}, {48, 105, 127, 248, 77, 241, 224, 247, 64, 156, 95, 182, 236, 170, 150, 162, 11, 205, 212, 94, 134, 133, 213, 110, 239, 250, 45, 35, 30, 26, 218, 99, 130, 69, 108, 143, 40, 211, 206, 132, 229, 7, 144, 2, 96}, {96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15}};

/**
 * @brief Computes 2 * PARAM_DELTA syndromes
 *
 * @param[out] syndromes Array of size 2 * PARAM_DELTA receiving the computed syndromes
 * @param[in] cdw Array of size PARAM_N1 storing the received vector
 */
static void compute_syndromes(uint16_t *syndromes, uint8_t *cdw) {
    for (size_t i = 0; i < 2 * PARAM_DELTA; ++i) {
        for (size_t j = 1; j < PARAM_N1; ++j) {
            syndromes[i] ^= PQCLEAN_HQC128_CLEAN_gf_mul(cdw[j], alpha_ij_pow[i][j - 1]);
        }
        syndromes[i] ^= cdw[0];
    }
}




/**
 * @brief Computes the error locator polynomial (ELP) sigma
 *
 * This is a constant time implementation of Berlekamp's simplified algorithm (see @cite lin1983error (Chapter 6 - BCH Codes). <br>
 * We use the letter p for rho which is initialized at -1. <br>
 * The array X_sigma_p represents the polynomial X^(mu-rho)*sigma_p(X). <br>
 * Instead of maintaining a list of sigmas, we update in place both sigma and X_sigma_p. <br>
 * sigma_copy serves as a temporary save of sigma in case X_sigma_p needs to be updated. <br>
 * We can properly correct only if the degree of sigma does not exceed PARAM_DELTA.
 * This means only the first PARAM_DELTA + 1 coefficients of sigma are of value
 * and we only need to save its first PARAM_DELTA - 1 coefficients.
 *
 * @returns the degree of the ELP sigma
 * @param[out] sigma Array of size (at least) PARAM_DELTA receiving the ELP
 * @param[in] syndromes Array of size (at least) 2*PARAM_DELTA storing the syndromes
 */
static uint16_t compute_elp(uint16_t *sigma, const uint16_t *syndromes) {
    uint16_t deg_sigma = 0;
    uint16_t deg_sigma_p = 0;
    uint16_t deg_sigma_copy = 0;
    uint16_t sigma_copy[PARAM_DELTA + 1] = {0};
    uint16_t X_sigma_p[PARAM_DELTA + 1] = {0, 1};
    uint16_t pp = (uint16_t) -1; // 2*rho
    uint16_t d_p = 1;
    uint16_t d = syndromes[0];

    uint16_t mask1, mask2, mask12;
    uint16_t deg_X, deg_X_sigma_p;
    uint16_t dd;
    uint16_t mu;

    uint16_t i;

    sigma[0] = 1;
    for (mu = 0; (mu < (2 * PARAM_DELTA)); ++mu) {
        // Save sigma in case we need it to update X_sigma_p
        memcpy(sigma_copy, sigma, 2 * (PARAM_DELTA));
        deg_sigma_copy = deg_sigma;

        dd = PQCLEAN_HQC128_CLEAN_gf_mul(d, PQCLEAN_HQC128_CLEAN_gf_inverse(d_p));

        for (i = 1; (i <= mu + 1) && (i <= PARAM_DELTA); ++i) {
            sigma[i] ^= PQCLEAN_HQC128_CLEAN_gf_mul(dd, X_sigma_p[i]);
        }

        deg_X = mu - pp;
        deg_X_sigma_p = deg_X + deg_sigma_p;

        // mask1 = 0xffff if(d != 0) and 0 otherwise
        mask1 = -((uint16_t) - d >> 15);

        // mask2 = 0xffff if(deg_X_sigma_p > deg_sigma) and 0 otherwise
        mask2 = -((uint16_t) (deg_sigma - deg_X_sigma_p) >> 15);

        // mask12 = 0xffff if the deg_sigma increased and 0 otherwise
        mask12 = mask1 & mask2;
        deg_sigma ^= mask12 & (deg_X_sigma_p ^ deg_sigma);

        if (mu == (2 * PARAM_DELTA - 1)) {
            break;
        }

        pp ^= mask12 & (mu ^ pp);
        d_p ^= mask12 & (d ^ d_p);
        for (i = PARAM_DELTA; i; --i) {
            X_sigma_p[i] = (mask12 & sigma_copy[i - 1]) ^ (~mask12 & X_sigma_p[i - 1]);
        }

        deg_sigma_p ^= mask12 & (deg_sigma_copy ^ deg_sigma_p);
        d = syndromes[mu + 1];

        for (i = 1; (i <= mu + 1) && (i <= PARAM_DELTA); ++i) {
            d ^= PQCLEAN_HQC128_CLEAN_gf_mul(sigma[i], syndromes[mu + 1 - i]);
        }
    }

    return deg_sigma;
}



/**
 * @brief Computes the polynomial z(x)
 *
 * See @cite lin1983error (Chapter 6 - BCH Codes) for more details.
 *
 * @param[out] z Array of PARAM_DELTA + 1 elements receiving the polynomial z(x)
 * @param[in] sigma Array of 2^PARAM_FFT elements storing the error locator polynomial
 * @param[in] degree Integer that is the degree of polynomial sigma
 * @param[in] syndromes Array of 2 * PARAM_DELTA storing the syndromes
 */
static void compute_z_poly(uint16_t *z, const uint16_t *sigma, uint16_t degree, const uint16_t *syndromes) {
    size_t i, j;
    uint16_t mask;

    z[0] = 1;

    for (i = 1; i < PARAM_DELTA + 1; ++i) {
        mask = -((uint16_t) (i - degree - 1) >> 15);
        z[i] = mask & sigma[i];
    }

    z[1] ^= syndromes[0];

    for (i = 2; i <= PARAM_DELTA; ++i) {
        mask = -((uint16_t) (i - degree - 1) >> 15);
        z[i] ^= mask & syndromes[i - 1];

        for (j = 1; j < i; ++j) {
            z[i] ^= mask & PQCLEAN_HQC128_CLEAN_gf_mul(sigma[j], syndromes[i - j - 1]);
        }
    }
}



/**
 * @brief Computes the error values
 *
 * See @cite lin1983error (Chapter 6 - BCH Codes) for more details.
 *
 * @param[out] error_values Array of PARAM_DELTA elements receiving the error values
 * @param[in] z Array of PARAM_DELTA + 1 elements storing the polynomial z(x)
 * @param[in] z_degree Integer that is the degree of polynomial z(x)
 * @param[in] error_compact Array of PARAM_DELTA + PARAM_N1 storing compact representation of the error
 */
static void compute_error_values(uint16_t *error_values, const uint16_t *z, const uint8_t *error) {
    uint16_t beta_j[PARAM_DELTA] = {0};
    uint16_t e_j[PARAM_DELTA] = {0};

    uint16_t delta_counter;
    uint16_t delta_real_value;
    uint16_t found;
    uint16_t mask1;
    uint16_t mask2;
    uint16_t tmp1;
    uint16_t tmp2;
    uint16_t inverse;
    uint16_t inverse_power_j;

    // Compute the beta_{j_i} page 31 of the documentation
    delta_counter = 0;
    for (size_t i = 0; i < PARAM_N1; i++) {
        found = 0;
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0
        for (size_t j = 0; j < PARAM_DELTA; j++) {
            mask2 = ~((uint16_t) (-((int32_t) j ^ delta_counter) >> 31)); // j == delta_counter
            beta_j[j] += mask1 & mask2 & gf_exp[i];
            found += mask1 & mask2 & 1;
        }
        delta_counter += found;
    }
    delta_real_value = delta_counter;

    // Compute the e_{j_i} page 31 of the documentation
    for (size_t i = 0; i < PARAM_DELTA; ++i) {
        tmp1 = 1;
        tmp2 = 1;
        inverse = PQCLEAN_HQC128_CLEAN_gf_inverse(beta_j[i]);
        inverse_power_j = 1;

        for (size_t j = 1; j <= PARAM_DELTA; ++j) {
            inverse_power_j = PQCLEAN_HQC128_CLEAN_gf_mul(inverse_power_j, inverse);
            tmp1 ^= PQCLEAN_HQC128_CLEAN_gf_mul(inverse_power_j, z[j]);
        }
        for (size_t k = 1; k < PARAM_DELTA; ++k) {
            tmp2 = PQCLEAN_HQC128_CLEAN_gf_mul(tmp2, (1 ^ PQCLEAN_HQC128_CLEAN_gf_mul(inverse, beta_j[(i + k) % PARAM_DELTA])));
        }
        mask1 = (uint16_t) (((int16_t) i - delta_real_value) >> 15); // i < delta_real_value
        e_j[i] = mask1 & PQCLEAN_HQC128_CLEAN_gf_mul(tmp1, PQCLEAN_HQC128_CLEAN_gf_inverse(tmp2));
    }

    // Place the delta e_{j_i} values at the right coordinates of the output vector
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
#endif

/*
static const uint32_t VandermoreMatrix_T[184] = {
    0x01010101, 0x01010101, 0x01010101, 0x01010101, 
    0xad478e01, 0x1b366cd8, 0xfae9cf83, 0x2c58b07d, 
    0x36d84701, 0x587de983, 0xf3eb8b16, 0x24907af5, 
    0xcf36ad01, 0xfb8b2c7d, 0x8a243df5, 0x59f2c356, 
    0x7d83d801, 0x90f5eb16, 0xb2ef5609, 0x640ee0a6, 
    0x2ce96c01, 0xac24f4eb, 0x1ca759ef, 0x9637ae64, 
    0x8b7d3601, 0xf25624f5, 0x576438a6, 0x91b3dbc4, 
    0xfb581b01, 0x53f2ac90, 0xab37820e, 0x55d5c6b3, 
    0xf5168301, 0x0ea6ef09, 0x7effc441, 0xa954aa73, 
    0x3d8bcf01, 0x82385956, 0xe491f1c4, 0x3b172129, 
    0x24ebe901, 0x3764a7ef, 0x4d7291ff, 0x1a7c33a9, 
    0x8af3fa01, 0xab571cb2, 0xda4de47e, 0xdf2281c5, 
    0x56f57d01, 0xb3c464a6, 0xc5a92973, 0x0f7f86ce, 
    0xc37ab001, 0xc6dbaee0, 0x813321aa, 0xb9cab186, 
    0xf2905801, 0xd5b3370e, 0x227c1754, 0xc1d2ca7f, 
    0x59242c01, 0x55919664, 0xdf1a3ba9, 0x60c1b90f, 
    0xa6091601, 0x5473ff41, 0xf0e2cecc, 0x26c023a1, 
    0xdd450b01, 0x4f92d7dc, 0x99d64493, 0x01984e0a, 
    0x38568b01, 0x172991c4, 0xba0fb6ce, 0x2c087535, 
    0x079bcb01, 0xec15e64b, 0x9fbe6b88, 0x24fa1d0c, 
    0x64efeb01, 0x7ca972ff, 0x9ca00fe2, 0x59f46c26, 
    0x82f2fb01, 0xd01755b3, 0x8fc1617f, 0x64c38b08, 
    0x57b2f301, 0x22c54d7e, 0x879cbaf0, 0x967012cf, 
    0x6ea2f701, 0x71ed8463, 0x040346bc, 0x91a5b2fb, 
    0xc4a6f501, 0x7fcea973, 0xcf2635a1, 0x55f1078a, 
    0x96a7f401, 0xfd1a2e72, 0xeb2060a0, 0xa9e63759, 
    0xdbe07a01, 0xca8633aa, 0x126c7523, 0x3b9a7b07, 
    0xf1383d01, 0x61b63b29, 0xf22ccd35, 0x1ab8736e, 
    0xb30e9001, 0xd27f7c54, 0x70f408c0, 0xdfc79af1, 
    0xfc8d4801, 0x14e76742, 0x41ac36c9, 0x0f0d6dbf, 
    0x91642401, 0xc10f1aa9, 0x96592c26, 0xb9df3b55, 
    0xbf191201, 0x942f115c, 0xe51cf580, 0xc11e6821, 
    0x73410901, 0xc0a1e2cc, 0xb7ae8a8e, 0x60de7166, 
    0xe4578a01, 0x8fbadfc5, 0x2996f2cf, 0x2646e73e, 
    0x92dc4501, 0x980ad693, 0x4fd7dd0b, 0x014e9944, 
    0x5537ac01, 0x74c1fd7c, 0x33e664f4, 0x2cb4a0df, 
    0x29c45601, 0x08350fce, 0x3e556e8a, 0x2440b578, 
    0xa8312b01, 0xd8275e68, 0x8884dbf9, 0x59833061, 
    0x154b9b01, 0xfa0cbe88, 0x5b2efc53, 0x64f34c50, 
    0x21dbc301, 0x8b75b986, 0xe73b7307, 0x965608b5, 
    0xa9ffef01, 0xf426a0e2, 0x5e6755ae, 0x91a7e960, 
    0xb8f6f901, 0x093a05a3, 0xd2111562, 0x5541f72d, 
    0x17b3f201, 0xc308c17f, 0x46dfb8f1, 0xa9db5640, 
    0x66e57901, 0xa2ad6abb, 0x4afdc53f, 0x3b635336, 
    0xc57eb201, 0x70cf9cf0, 0x065e3eb7, 0x1aa4198b, 
    0x3b915901, 0x642c600f, 0x26b91a55, 0xdfa99624
};
*/

static const uint32_t VandermoreMatrix_rdx16[192] = {
    0x0000000f, 0xc280ce6b, 0x268a2cc5, 0xc862cead, 0x6c8cc82b, 0x26e0ec03, 0x8266a6cb, 0x8a8ee0cd, 0xc88a0a6d, 0x6488ec6f, 0x66e06827, 0xe666c0e3, 0x2e6e4e87, 0xca6620a9, 0xac8e4081, 0x0868c609, 
    0x0000000f, 0x156bd6a8, 0x5e6ae217, 0xa9d9f5ab, 0xe66b4356, 0xd5f49e11, 0x9bdb07a1, 0x7a6b44a8, 0x6460dee5, 0x85276194, 0x53f81fdd, 0x9297c6ba, 0xb6d80798, 0xfb526ce2, 0xa86e1d7e, 0x691e0d4a, 
    0x0000000f, 0xfea8e4d7, 0xeaa961fe, 0x9175c785, 0xa6ac36e3, 0x353cdb37, 0x1a7c4b98, 0xc0a293fa, 0x6faaaea3, 0xfdc60906, 0x5b37993d, 0xcf59b792, 0xa971681f, 0x9072103e, 0x03a727c4, 0x985ef8ed, 
    0x0000000f, 0x43d7b901, 0x35d74941, 0x6f2b817c, 0x5cd02b30, 0xc17a4fe2, 0xf52e416e, 0x7e1b0e5b, 0xc3d3c259, 0x00da92af, 0x127eb2c5, 0xd92889ad, 0x5322bffa, 0xf2cca13c, 0xef135179, 0x6ac45117, 
    0x0000000f, 0xa801e1fe, 0x860a96ad, 0x54366d9c, 0x6f04c689, 0xfd971dd7, 0x413d285b, 0x8d83cf66, 0xf304b364, 0x3da297b8, 0xd59caef2, 0xf13fe8ac, 0x1032ca43, 0x063ace5b, 0xd98ce080, 0x00b02b14, 
    0x0000000f, 0xfffe58ff, 0xf5f611f1, 0x4b928be2, 0x5bf733fd, 0x87152e16, 0xba949e4f, 0x4dbd1128, 0xb7fecc54, 0x2c6239b5, 0x7412228c, 0x70015f66, 0xaa98e3ba, 0xc19a51e4, 0xd7b3094c, 0xa54b1629, 
    0x0000000f, 0x57ffad29, 0x72f88554, 0xfb8e256a, 0x25f5a079, 0x10fad4eb, 0xb18c5bf0, 0x6fd5753d, 0x50f8082b, 0x8e738b0b, 0x02f0ab10, 0xddca51f5, 0x1c84dab6, 0x0aee25ff, 0xf5dba461, 0x416be1a7, 
    0x0000000f, 0x14294280, 0x4228d71c, 0xf4dfcc9b, 0x2d29a942, 0x055b356c, 0x45d1bfff, 0xb30d0c12, 0xd5204224, 0x88ec3e48, 0x535bed02, 0xad930a2c, 0x53d7464f, 0x1ce3f436, 0x330aefb7, 0xb923755f, 
    0x0000000f, 0xc280da43, 0x278e2ecd, 0x6d5b0c39, 0x7d8a9d26, 0xff9ece7c, 0xd8546069, 0x798a6c44, 0xdc887d76, 0x7ec63a52, 0xfa961ff2, 0x12d963bf, 0x8b5ccfd8, 0xe5d32447, 0x97893f70, 0x01d08c42, 
    0x0000000f, 0xa94361fe, 0x9b4ca1ae, 0x40d2520b, 0xb5407990, 0xa4aef779, 0x08d92f4d, 0xe626f0a8, 0x564658be, 0x6fbe73dc, 0x41a9a7a9, 0xca9379ba, 0x87db180c, 0x1255a301, 0x6b2a08e6, 0x863c4b84, 
    0x0000000f, 0x7ffed13d, 0xfdf71d7b, 0x58d65aa3, 0xd3faa2fe, 0x3a908486, 0x80d4e958, 0xb39141a6, 0x36fa99da, 0xa9272d0d, 0xab9eab3f, 0x4f0593f0, 0x03d84c84, 0x5f5e16ed, 0x3e9221bf, 0x34ada48f, 
    0x00000005, 0x14154014, 0x41151114, 0x11444145, 0x11114444, 0x14411100, 0x11405511, 0x05510400, 0x14150515, 0x00155150, 0x41455114, 0x00400550, 0x10401411, 0x40051441, 0x50504104, 0x44055454,
};

#if defined(_M4_ASM_)
void radix_16_matrix_vector_46x16(uint32_t *in, const uint32_t *mat, uint32_t *vec);
#else
void radix_16_matrix_vector_46x16(uint32_t *in, const uint32_t *mat, uint32_t *vec) {
    for(int i=0;i<(11);i++) { in[i] = radix_16_matrix_vector(&mat[16*i], vec, 16); } 
    in[11] = radix_16_matrix_vector_half(&mat[16*11], vec, 16);
}
#endif

static void compute_roots_new(uint8_t *error, uint8_t *sigma) {
    
    /*
    uint8_t * vande_u8 = (uint8_t *)VandermoreMatrix_T;
    for(int i = 0; i < PARAM_N1; i++) {
        uint8_t value = 0;
        for(int j = 0; j < PARAM_DELTA + 1; j++)
        {
            value ^= PQCLEAN_HQC128_CLEAN_gf_mul(sigma[j], vande_u8[i * (PARAM_DELTA + 1) + j]);
        }
        error[i] = (value == 0);
    }
    */

    uint32_t _sigma_rdx16[4];
    memcpy( _sigma_rdx16 , sigma , 16 );
    for(int i=0;i<4;i++) { _sigma_rdx16[i] = rdx16_from_bitseq( _sigma_rdx16[i] ); }
    uint32_t sigma_rdx16[16];
    for(int i=0; i<4; i++) {
        sigma_rdx16[4*i] = _sigma_rdx16[i]&_NORM_RDX16_;
        sigma_rdx16[4*i+1] = (_sigma_rdx16[i]>>2)&_NORM_RDX16_;
        sigma_rdx16[4*i+2] = (_sigma_rdx16[i]>>1)&_NORM_RDX16_;
        sigma_rdx16[4*i+3] = (_sigma_rdx16[i]>>3)&_NORM_RDX16_;
    }
        
    uint32_t product[12] = {0};
    radix_16_matrix_vector_46x16(&product[0] ,&VandermoreMatrix_rdx16[0], sigma_rdx16);
    for(int i=0; i<11; i++) { 
        error[4*i] = ((product[i]&0x11111111) == 0);
        error[4*i+1] = ((product[i]&0x44444444) == 0);
        error[4*i+2] = ((product[i]&0x22222222) == 0);
        error[4*i+3] = ((product[i]&0x88888888) == 0);
    }
    error[44] = ((product[11]&0x11111111) == 0);
    error[45] = ((product[11]&0x44444444) == 0);
    
    /*
    for(int i=0; i<11; i++) { // i in range(11)
        uint32_t r = radix_16_matrix_vector(&VandermoreMatrix_rdx16[16*i], sigma_rdx16);
        error[4*i] = ((r&0x11111111) == 0);
        error[4*i+1] = ((r&0x44444444) == 0);
        error[4*i+2] = ((r&0x22222222) == 0);
        error[4*i+3] = ((r&0x88888888) == 0);
    }
    { // i = 11
        uint32_t r = radix_16_matrix_vector(&VandermoreMatrix_rdx16[176], sigma_rdx16);
        error[44] = ((r&0x11111111) == 0);
        error[45] = ((r&0x44444444) == 0);
    }
    */
}
/**
 * @brief Computes the error polynomial error from the error locator polynomial sigma
 *
 * See function PQCLEAN_HQC128_CLEAN_fft for more details.
 *
 * @param[out] error Array of 2^PARAM_M elements receiving the error polynomial
 * @param[out] error_compact Array of PARAM_DELTA + PARAM_N1 elements receiving a compact representation of the vector error
 * @param[in] sigma Array of 2^PARAM_FFT elements storing the error locator polynomial
 */
static void compute_roots(uint8_t *error, uint16_t *sigma) {
    uint16_t w[1 << PARAM_M] = {0};

    PQCLEAN_HQC128_CLEAN_fft(w, sigma, PARAM_DELTA + 1);
    PQCLEAN_HQC128_CLEAN_fft_retrieve_error_poly(error, w);
}

#if 1
static void correct_errors(uint8_t *cdw, const uint8_t *error_values) {
    for (size_t i = 0; i < PARAM_N1; ++i) {
        cdw[i] ^= error_values[i];
    }
}
#else

/**
 * @brief Correct the errors
 *
 * @param[out] cdw Array of PARAM_N1 elements receiving the corrected vector
 * @param[in] error Array of the error vector
 * @param[in] error_values Array of PARAM_DELTA elements storing the error values
 */
static void correct_errors(uint8_t *cdw, const uint16_t *error_values) {
    for (size_t i = 0; i < PARAM_N1; ++i) {
        cdw[i] ^= error_values[i];
    }
}
#endif



/**
 * @brief Decodes the received word
 *
 * This function relies on six steps:
 *    <ol>
 *    <li> The first step, is the computation of the 2*PARAM_DELTA syndromes.
 *    <li> The second step is the computation of the error-locator polynomial sigma.
 *    <li> The third step, done by additive FFT, is finding the error-locator numbers by calculating the roots of the polynomial sigma and takings their inverses.
 *    <li> The fourth step, is the polynomial z(x).
 *    <li> The fifth step, is the computation of the error values.
 *    <li> The sixth step is the correction of the errors in the received polynomial.
 *    </ol>
 * For a more complete picture on Reed-Solomon decoding, see Shu. Lin and Daniel J. Costello in Error Control Coding: Fundamentals and Applications @cite lin1983error
 *
 * @param[out] msg Array of size VEC_K_SIZE_64 receiving the decoded message
 * @param[in] cdw Array of size VEC_N1_SIZE_64 storing the received word
 */
void PQCLEAN_HQC128_CLEAN_reed_solomon_decode(uint8_t *msg, uint8_t *cdw) {
#if 1

    uint8_t syndromes_new[2 * PARAM_DELTA] = {0};
    uint8_t sigma_new[1 << PARAM_FFT] = {0};
    uint16_t sigma_new_16[1 << PARAM_FFT] = {0};
    uint8_t omega[1 << PARAM_FFT] = {0};
    uint8_t error_values_new[PARAM_N1] = {0};
    uint8_t error_new[1 << PARAM_M] = {0};
    uint8_t error[1 << PARAM_M] = {0};
    uint8_t mask = 0;

    // Calculate the 2*PARAM_DELTA syndromes
    compute_syndromes(syndromes_new, cdw);

    // Compute the error locator polynomial sigma
    // Sigma's degree is at most PARAM_DELTA but the FFT requires the extra room 
    compute_elp_divstep(omega, sigma_new, syndromes_new);
    for(int i = 0; i < 16; i++)
    {
        sigma_new_16[i] = sigma_new[i];
        mask ^= sigma_new[i];
    }

    // Compute the error polynomial error
    compute_roots_new(error_new, sigma_new);
    // compute_roots(error_new, sigma_new_16);
    // error_new[0] = (mask == 0);
    // for(int i = 0; i < (1 << PARAM_M); i++)
    // {
    //     if(error[i] != error_new[i])
    //     {
    //         char buffer[100];
    //         sprintf(buffer, "index %d: error=%d, error_new=%d\n", i, error[i], error_new[i]);
    //         hal_send_str(buffer);
    //     }
    // }


    // Compute the error values
    compute_error_values_new(error_values_new, omega, sigma_new, error_new);
    correct_errors(cdw, error_values_new);
#else
    uint16_t syndromes[2 * PARAM_DELTA] = {0};
    uint16_t sigma[1 << PARAM_FFT] = {0};
    uint8_t error[1 << PARAM_M] = {0};
    uint16_t z[PARAM_N1] = {0};
    uint16_t error_values[PARAM_N1] = {0};
    uint16_t deg;

    // Calculate the 2*PARAM_DELTA syndromes
    compute_syndromes(syndromes, cdw);

    // Compute the error locator polynomial sigma
    // Sigma's degree is at most PARAM_DELTA but the FFT requires the extra room
    deg = compute_elp(sigma, syndromes);

    // Compute the error polynomial error
    compute_roots(error, sigma);

    // Compute the polynomial z(x)
    compute_z_poly(z, sigma, deg, syndromes);

    // Compute the error values
    compute_error_values(error_values, z, error);

    // Correct the errors
    correct_errors(cdw, error_values);
#endif

    // Retrieve the message from the decoded codeword
    memcpy(msg, cdw + (PARAM_G - 1), PARAM_K);

}
