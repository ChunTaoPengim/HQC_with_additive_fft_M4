#include "fft.h"
#include "gf.h"
#include "parameters.h"
#include "reed_solomon.h"
#include "gf256_rdx16.h"
#include "gf256.h"
#include <stdint.h>
#include <string.h>
// #include "hal.h"
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
void PQCLEAN_HQC192_CLEAN_reed_solomon_encode(uint8_t *cdw, const uint8_t *msg) {
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
void PQCLEAN_HQC192_CLEAN_reed_solomon_encode(uint8_t *cdw, const uint8_t *msg) {
    uint8_t gate_value = 0;

    uint16_t tmp[PARAM_G] = {0};
    uint16_t PARAM_RS_POLY [] = {RS_POLY_COEFS};

    memset(cdw, 0, PARAM_N1);

    for (size_t i = 0; i < PARAM_K; ++i) {
        gate_value = msg[PARAM_K - 1 - i] ^ cdw[PARAM_N1 - PARAM_K - 1];

        for (size_t j = 0; j < PARAM_G; ++j) {
            tmp[j] = PQCLEAN_HQC192_CLEAN_gf_mul(gate_value, PARAM_RS_POLY[j]);
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

#elif 1

static const uint8_t alpha_ij_pow_trans[55][32] ={
    {2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157},
    {4, 16, 64, 29, 116, 205, 19, 76, 45, 180, 234, 143, 6, 24, 96, 157, 78, 37, 148, 106, 181, 238, 159, 70, 5, 20, 80, 93, 105, 185, 222, 95},
    {8, 64, 58, 205, 38, 45, 117, 143, 12, 96, 39, 37, 53, 181, 193, 70, 10, 80, 186, 185, 161, 97, 47, 101, 15, 120, 231, 107, 127, 223, 182, 217},
    {16, 29, 205, 76, 180, 143, 24, 157, 37, 106, 238, 70, 20, 93, 185, 95, 153, 101, 30, 253, 107, 254, 91, 217, 17, 13, 208, 129, 248, 59, 151, 133},
    {32, 116, 38, 180, 3, 96, 156, 106, 193, 5, 160, 185, 190, 94, 15, 253, 214, 223, 226, 17, 26, 103, 124, 59, 51, 46, 169, 132, 77, 85, 114, 230},
    {64, 205, 45, 143, 96, 37, 181, 70, 80, 185, 97, 101, 120, 107, 223, 217, 68, 208, 62, 59, 102, 184, 33, 168, 85, 228, 191, 252, 241, 150, 110, 130},
    {128, 19, 117, 24, 156, 181, 140, 93, 161, 94, 60, 107, 163, 67, 26, 129, 147, 102, 109, 132, 41, 57, 209, 252, 255, 98, 87, 200, 224, 89, 155, 18},
    {29, 76, 143, 157, 106, 70, 93, 95, 101, 253, 254, 217, 13, 129, 59, 133, 79, 168, 73, 230, 252, 227, 149, 130, 28, 81, 195, 18, 247, 44, 27, 2},
    {58, 45, 12, 37, 193, 80, 161, 101, 231, 223, 134, 208, 237, 102, 169, 168, 146, 191, 179, 150, 87, 7, 166, 195, 36, 251, 125, 173, 64, 38, 143, 39},
    {116, 180, 96, 106, 5, 185, 94, 253, 223, 17, 103, 59, 46, 132, 85, 230, 215, 150, 174, 28, 89, 172, 244, 44, 108, 32, 38, 3, 156, 193, 160, 190},
    {232, 234, 39, 238, 160, 97, 60, 254, 134, 103, 118, 184, 84, 57, 145, 227, 220, 7, 162, 172, 245, 176, 71, 58, 180, 192, 181, 40, 95, 15, 177, 175},
    {205, 143, 37, 70, 185, 101, 107, 217, 208, 59, 184, 168, 228, 252, 150, 130, 221, 195, 61, 44, 173, 58, 117, 39, 193, 186, 47, 231, 182, 26, 237, 23},
    {135, 6, 53, 20, 190, 120, 163, 13, 237, 46, 84, 228, 229, 98, 100, 81, 69, 251, 131, 32, 45, 192, 238, 186, 94, 187, 217, 189, 236, 169, 82, 209},
    {19, 24, 181, 93, 94, 107, 67, 129, 102, 132, 57, 252, 98, 200, 89, 18, 11, 173, 232, 3, 53, 40, 194, 231, 226, 189, 197, 158, 170, 145, 75, 25},
    {38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36},
    {76, 157, 70, 95, 253, 217, 129, 133, 168, 230, 227, 130, 81, 18, 44, 2, 152, 39, 140, 190, 231, 175, 31, 23, 77, 209, 219, 25, 162, 36, 88, 4},
    {152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78},
    {45, 37, 80, 101, 223, 208, 102, 168, 191, 150, 7, 195, 251, 173, 38, 39, 10, 47, 127, 26, 197, 21, 115, 219, 100, 242, 245, 54, 205, 96, 70, 97},
    {90, 148, 186, 30, 226, 62, 109, 73, 179, 174, 162, 61, 131, 232, 96, 140, 153, 127, 52, 51, 168, 99, 98, 56, 172, 22, 8, 234, 212, 185, 240, 67},
    {180, 106, 185, 253, 17, 59, 132, 230, 150, 28, 172, 44, 32, 3, 193, 190, 214, 26, 51, 77, 145, 55, 167, 36, 233, 116, 96, 5, 94, 223, 103, 46},
    {117, 181, 161, 107, 26, 102, 41, 252, 87, 89, 245, 173, 45, 53, 185, 231, 68, 197, 168, 145, 110, 166, 61, 54, 38, 37, 186, 120, 134, 59, 21, 191},
    {234, 238, 97, 254, 103, 184, 57, 227, 7, 172, 176, 58, 192, 40, 15, 175, 147, 21, 99, 55, 166, 122, 216, 45, 106, 222, 107, 52, 133, 85, 123, 50},
    {201, 159, 47, 91, 124, 33, 209, 149, 166, 244, 71, 117, 238, 194, 223, 31, 79, 115, 98, 167, 61, 216, 90, 181, 190, 254, 206, 218, 213, 150, 224, 72},
    {143, 70, 101, 217, 59, 168, 252, 130, 195, 44, 58, 39, 186, 231, 26, 23, 146, 219, 56, 36, 54, 45, 181, 97, 223, 62, 33, 191, 110, 89, 251, 8},
    {3, 5, 15, 17, 51, 85, 255, 28, 36, 108, 180, 193, 94, 226, 59, 77, 215, 100, 172, 233, 38, 106, 190, 223, 124, 132, 145, 174, 239, 44, 116, 156},
    {6, 20, 120, 13, 46, 228, 98, 81, 251, 32, 192, 186, 187, 189, 169, 209, 220, 242, 22, 116, 37, 222, 254, 62, 132, 63, 130, 43, 250, 38, 212, 194},
    {12, 80, 231, 208, 169, 191, 87, 195, 125, 38, 181, 47, 217, 197, 85, 219, 221, 245, 8, 96, 186, 107, 206, 33, 145, 130, 86, 207, 45, 193, 101, 134},
    {24, 93, 107, 129, 132, 252, 200, 18, 173, 3, 40, 231, 189, 158, 145, 25, 69, 54, 234, 5, 120, 52, 218, 191, 174, 43, 207, 90, 35, 15, 136, 92},
    {48, 105, 127, 248, 77, 241, 224, 247, 64, 156, 95, 182, 236, 170, 150, 162, 11, 205, 212, 94, 134, 133, 213, 110, 239, 250, 45, 35, 30, 26, 218, 99},
    {96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100},
    {192, 222, 182, 151, 114, 110, 155, 27, 143, 160, 177, 237, 82, 75, 89, 88, 152, 70, 240, 103, 21, 123, 224, 251, 116, 212, 101, 136, 218, 145, 200, 144},
    {157, 95, 217, 133, 230, 130, 18, 2, 39, 190, 175, 23, 209, 25, 36, 4, 78, 97, 67, 46, 191, 50, 72, 8, 156, 194, 134, 92, 99, 100, 144, 16},
    {39, 97, 134, 184, 145, 7, 245, 58, 181, 15, 208, 21, 241, 166, 44, 45, 10, 107, 237, 85, 196, 195, 54, 12, 185, 182, 102, 115, 130, 36, 8, 37},
    {78, 153, 68, 79, 215, 221, 11, 152, 10, 214, 147, 146, 220, 69, 1, 78, 153, 68, 79, 215, 221, 11, 152, 10, 214, 147, 146, 220, 69, 1, 78, 153},
    {156, 94, 26, 132, 255, 89, 233, 3, 185, 226, 46, 145, 28, 235, 38, 5, 214, 59, 114, 174, 36, 32, 106, 15, 103, 77, 150, 239, 108, 96, 190, 17},
    {37, 101, 208, 168, 150, 195, 173, 39, 47, 26, 21, 219, 242, 54, 96, 97, 68, 33, 241, 89, 207, 12, 161, 134, 169, 179, 166, 125, 143, 185, 217, 184},
    {74, 137, 206, 82, 55, 138, 16, 212, 120, 124, 73, 87, 72, 29, 193, 211, 147, 228, 25, 244, 205, 140, 177, 197, 230, 141, 251, 76, 40, 223, 204, 198},
    {148, 30, 62, 73, 174, 61, 232, 140, 127, 51, 99, 56, 22, 234, 185, 67, 79, 241, 121, 108, 39, 188, 189, 41, 55, 9, 64, 238, 211, 59, 183, 200},
    {53, 120, 237, 228, 100, 251, 45, 186, 217, 169, 241, 242, 173, 37, 15, 62, 146, 130, 245, 38, 80, 182, 184, 179, 89, 54, 39, 101, 206, 85, 87, 61},
    {106, 253, 59, 230, 28, 44, 3, 190, 26, 77, 55, 36, 116, 5, 223, 46, 215, 89, 108, 156, 15, 124, 114, 100, 235, 180, 185, 17, 132, 150, 172, 32},
    {212, 211, 197, 198, 167, 207, 157, 202, 62, 114, 200, 139, 201, 95, 26, 154, 220, 61, 19, 160, 217, 158, 171, 86, 32, 159, 127, 133, 229, 89, 216, 74},
    {181, 107, 102, 252, 89, 173, 53, 231, 197, 145, 166, 54, 37, 120, 59, 191, 221, 207, 39, 15, 237, 115, 56, 125, 96, 101, 62, 228, 7, 44, 12, 47},
    {119, 177, 23, 123, 239, 8, 159, 225, 184, 255, 43, 64, 140, 91, 169, 171, 69, 58, 20, 226, 33, 49, 18, 205, 160, 67, 21, 149, 144, 38, 105, 34},
    {238, 254, 184, 227, 172, 58, 40, 175, 21, 55, 122, 45, 222, 52, 85, 50, 11, 12, 188, 124, 115, 224, 131, 37, 253, 151, 252, 121, 2, 193, 225, 109},
    {193, 223, 169, 150, 36, 38, 185, 26, 85, 100, 44, 96, 15, 59, 145, 89, 1, 193, 223, 169, 150, 36, 38, 185, 26, 85, 100, 44, 96, 15, 59, 145},
    {159, 91, 33, 149, 244, 117, 194, 31, 115, 167, 216, 181, 254, 218, 150, 72, 152, 161, 189, 114, 56, 131, 148, 107, 46, 227, 138, 135, 210, 26, 170, 141},
    {35, 113, 21, 165, 235, 12, 137, 118, 252, 239, 128, 80, 34, 82, 100, 176, 78, 231, 133, 255, 138, 19, 111, 208, 114, 112, 54, 212, 254, 169, 98, 122},
    {70, 217, 168, 130, 44, 39, 231, 23, 219, 36, 45, 97, 62, 191, 89, 8, 10, 134, 41, 100, 125, 37, 107, 184, 150, 61, 117, 47, 237, 145, 242, 64},
    {140, 67, 41, 200, 233, 53, 254, 158, 110, 235, 48, 120, 204, 227, 36, 90, 153, 237, 63, 239, 58, 105, 104, 228, 167, 142, 70, 175, 154, 100, 250, 148},
    {5, 17, 85, 28, 108, 193, 226, 77, 100, 233, 106, 223, 132, 174, 44, 156, 214, 169, 55, 235, 96, 253, 46, 150, 244, 3, 15, 51, 255, 36, 180, 94},
    {10, 68, 146, 221, 1, 10, 68, 146, 221, 1, 10, 68, 146, 221, 1, 10, 68, 146, 221, 1, 10, 68, 146, 221, 1, 10, 68, 146, 221, 1, 10, 68},
    {20, 13, 228, 81, 32, 186, 189, 209, 242, 116, 222, 62, 63, 43, 38, 194, 147, 179, 9, 180, 101, 151, 227, 61, 3, 60, 23, 49, 243, 96, 211, 218},
    {40, 52, 115, 121, 116, 161, 248, 229, 138, 180, 202, 102, 75, 247, 96, 187, 79, 87, 176, 106, 182, 154, 14, 173, 5, 136, 228, 162, 128, 185, 31, 63},
    {80, 208, 191, 195, 38, 47, 197, 219, 245, 96, 107, 33, 130, 207, 193, 134, 146, 166, 64, 185, 62, 252, 138, 117, 15, 23, 196, 139, 37, 223, 168, 7},
    {160, 103, 145, 172, 180, 15, 46, 55, 44, 106, 226, 85, 167, 32, 185, 124, 215, 36, 3, 253, 169, 174, 233, 193, 17, 114, 89, 116, 190, 59, 255, 244}
};

/**
 * @brief Computes 2 * PARAM_DELTA syndromes
 *
 * @param[out] syndromes Array of size 2 * PARAM_DELTA receiving the computed syndromes
 * @param[in] cdw Array of size PARAM_N1 storing the received vector
 */
static void compute_syndromes(uint8_t *syndromes, uint8_t *cdw) {

    memset(syndromes, cdw[0], 2 * PARAM_DELTA);
    for (size_t j = 1; j < PARAM_N1; ++j) {
        gf256_madd((uint32_t *)syndromes, (uint32_t*)alpha_ij_pow_trans[j-1], (uint32_t)cdw[j]);
    }
}
#else

static const uint16_t alpha_ij_pow [32][55] = {{2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35, 70, 140, 5, 10, 20, 40, 80, 160}, {4, 16, 64, 29, 116, 205, 19, 76, 45, 180, 234, 143, 6, 24, 96, 157, 78, 37, 148, 106, 181, 238, 159, 70, 5, 20, 80, 93, 105, 185, 222, 95, 97, 153, 94, 101, 137, 30, 120, 253, 211, 107, 177, 254, 223, 91, 113, 217, 67, 17, 68, 13, 52, 208, 103}, {8, 64, 58, 205, 38, 45, 117, 143, 12, 96, 39, 37, 53, 181, 193, 70, 10, 80, 186, 185, 161, 97, 47, 101, 15, 120, 231, 107, 127, 223, 182, 217, 134, 68, 26, 208, 206, 62, 237, 59, 197, 102, 23, 184, 169, 33, 21, 168, 41, 85, 146, 228, 115, 191, 145}, {16, 29, 205, 76, 180, 143, 24, 157, 37, 106, 238, 70, 20, 93, 185, 95, 153, 101, 30, 253, 107, 254, 91, 217, 17, 13, 208, 129, 248, 59, 151, 133, 184, 79, 132, 168, 82, 73, 228, 230, 198, 252, 123, 227, 150, 149, 165, 130, 200, 28, 221, 81, 121, 195, 172}, {32, 116, 38, 180, 3, 96, 156, 106, 193, 5, 160, 185, 190, 94, 15, 253, 214, 223, 226, 17, 26, 103, 124, 59, 51, 46, 169, 132, 77, 85, 114, 230, 145, 215, 255, 150, 55, 174, 100, 28, 167, 89, 239, 172, 36, 244, 235, 44, 233, 108, 1, 32, 116, 38, 180}, {64, 205, 45, 143, 96, 37, 181, 70, 80, 185, 97, 101, 120, 107, 223, 217, 68, 208, 62, 59, 102, 184, 33, 168, 85, 228, 191, 252, 241, 150, 110, 130, 7, 221, 89, 195, 138, 61, 251, 44, 207, 173, 8, 58, 38, 117, 12, 39, 53, 193, 10, 186, 161, 47, 15}, {128, 19, 117, 24, 156, 181, 140, 93, 161, 94, 60, 107, 163, 67, 26, 129, 147, 102, 109, 132, 41, 57, 209, 252, 255, 98, 87, 200, 224, 89, 155, 18, 245, 11, 233, 173, 16, 232, 45, 3, 157, 53, 159, 40, 185, 194, 137, 231, 254, 226, 68, 189, 248, 197, 46}, {29, 76, 143, 157, 106, 70, 93, 95, 101, 253, 254, 217, 13, 129, 59, 133, 79, 168, 73, 230, 252, 227, 149, 130, 28, 81, 195, 18, 247, 44, 27, 2, 58, 152, 3, 39, 212, 140, 186, 190, 202, 231, 225, 175, 26, 31, 118, 23, 158, 77, 146, 209, 229, 219, 55}, {58, 45, 12, 37, 193, 80, 161, 101, 231, 223, 134, 208, 237, 102, 169, 168, 146, 191, 179, 150, 87, 7, 166, 195, 36, 251, 125, 173, 64, 38, 143, 39, 181, 10, 185, 47, 120, 127, 217, 26, 62, 197, 184, 21, 85, 115, 252, 219, 110, 100, 221, 242, 138, 245, 44}, {116, 180, 96, 106, 5, 185, 94, 253, 223, 17, 103, 59, 46, 132, 85, 230, 215, 150, 174, 28, 89, 172, 244, 44, 108, 32, 38, 3, 156, 193, 160, 190, 15, 214, 226, 26, 124, 51, 169, 77, 114, 145, 255, 55, 100, 167, 239, 36, 235, 233, 1, 116, 180, 96, 106}, {232, 234, 39, 238, 160, 97, 60, 254, 134, 103, 118, 184, 84, 57, 145, 227, 220, 7, 162, 172, 245, 176, 71, 58, 180, 192, 181, 40, 95, 15, 177, 175, 208, 147, 46, 21, 73, 99, 241, 55, 200, 166, 43, 122, 44, 216, 128, 45, 48, 106, 10, 222, 202, 107, 226}, {205, 143, 37, 70, 185, 101, 107, 217, 208, 59, 184, 168, 228, 252, 150, 130, 221, 195, 61, 44, 173, 58, 117, 39, 193, 186, 47, 231, 182, 26, 237, 23, 21, 146, 145, 219, 87, 56, 242, 36, 139, 54, 64, 45, 96, 181, 80, 97, 120, 223, 68, 62, 102, 33, 85}, {135, 6, 53, 20, 190, 120, 163, 13, 237, 46, 84, 228, 229, 98, 100, 81, 69, 251, 131, 32, 45, 192, 238, 186, 94, 187, 217, 189, 236, 169, 82, 209, 241, 220, 28, 242, 72, 22, 173, 116, 201, 37, 140, 222, 15, 254, 34, 62, 204, 132, 146, 63, 75, 130, 167}, {19, 24, 181, 93, 94, 107, 67, 129, 102, 132, 57, 252, 98, 200, 89, 18, 11, 173, 232, 3, 53, 40, 194, 231, 226, 189, 197, 158, 170, 145, 75, 25, 166, 69, 235, 54, 29, 234, 37, 5, 95, 120, 91, 52, 59, 218, 82, 191, 227, 174, 221, 43, 247, 207, 32}, {38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185, 15, 223, 26, 59, 169, 85, 145, 150, 100, 89, 36, 44, 1, 38, 96, 193, 185}, {76, 157, 70, 95, 253, 217, 129, 133, 168, 230, 227, 130, 81, 18, 44, 2, 152, 39, 140, 190, 231, 175, 31, 23, 77, 209, 219, 25, 162, 36, 88, 4, 45, 78, 5, 97, 211, 67, 62, 46, 154, 191, 171, 50, 89, 72, 176, 8, 90, 156, 10, 194, 187, 134, 124}, {152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215, 220, 221, 69, 11, 1, 152, 78, 10, 153, 214, 68, 147, 79, 146, 215}, {45, 37, 80, 101, 223, 208, 102, 168, 191, 150, 7, 195, 251, 173, 38, 39, 10, 47, 127, 26, 197, 21, 115, 219, 100, 242, 245, 54, 205, 96, 70, 97, 107, 68, 59, 33, 228, 241, 130, 89, 61, 207, 58, 12, 193, 161, 231, 134, 237, 169, 146, 179, 87, 166, 36}, {90, 148, 186, 30, 226, 62, 109, 73, 179, 174, 162, 61, 131, 232, 96, 140, 153, 127, 52, 51, 168, 99, 98, 56, 172, 22, 8, 234, 212, 185, 240, 67, 237, 79, 114, 241, 25, 121, 245, 108, 19, 39, 20, 188, 223, 189, 133, 41, 63, 55, 221, 9, 176, 64, 3}, {180, 106, 185, 253, 17, 59, 132, 230, 150, 28, 172, 44, 32, 3, 193, 190, 214, 26, 51, 77, 145, 55, 167, 36, 233, 116, 96, 5, 94, 223, 103, 46, 85, 215, 174, 89, 244, 108, 38, 156, 160, 15, 226, 124, 169, 114, 255, 100, 239, 235, 1, 180, 106, 185, 253}, {117, 181, 161, 107, 26, 102, 41, 252, 87, 89, 245, 173, 45, 53, 185, 231, 68, 197, 168, 145, 110, 166, 61, 54, 38, 37, 186, 120, 134, 59, 21, 191, 196, 221, 36, 207, 205, 39, 80, 15, 217, 237, 33, 115, 150, 56, 138, 125, 58, 96, 10, 101, 182, 62, 169}, {234, 238, 97, 254, 103, 184, 57, 227, 7, 172, 176, 58, 192, 40, 15, 175, 147, 21, 99, 55, 166, 122, 216, 45, 106, 222, 107, 52, 133, 85, 123, 50, 195, 11, 32, 12, 140, 188, 182, 124, 158, 115, 49, 224, 36, 131, 19, 37, 105, 253, 68, 151, 154, 252, 174}, {201, 159, 47, 91, 124, 33, 209, 149, 166, 244, 71, 117, 238, 194, 223, 31, 79, 115, 98, 167, 61, 216, 90, 181, 190, 254, 206, 218, 213, 150, 224, 72, 54, 152, 106, 161, 177, 189, 184, 114, 171, 56, 18, 131, 38, 148, 111, 107, 104, 46, 146, 227, 14, 138, 233}, {143, 70, 101, 217, 59, 168, 252, 130, 195, 44, 58, 39, 186, 231, 26, 23, 146, 219, 56, 36, 54, 45, 181, 97, 223, 62, 33, 191, 110, 89, 251, 8, 12, 10, 15, 134, 197, 41, 179, 100, 86, 125, 205, 37, 185, 107, 208, 184, 228, 150, 221, 61, 173, 117, 193}, {3, 5, 15, 17, 51, 85, 255, 28, 36, 108, 180, 193, 94, 226, 59, 77, 215, 100, 172, 233, 38, 106, 190, 223, 124, 132, 145, 174, 239, 44, 116, 156, 185, 214, 103, 169, 230, 55, 89, 235, 32, 96, 160, 253, 26, 46, 114, 150, 167, 244, 1, 3, 5, 15, 17}, {6, 20, 120, 13, 46, 228, 98, 81, 251, 32, 192, 186, 187, 189, 169, 209, 220, 242, 22, 116, 37, 222, 254, 62, 132, 63, 130, 43, 250, 38, 212, 194, 182, 147, 77, 179, 141, 9, 54, 180, 159, 101, 67, 151, 85, 227, 112, 61, 142, 3, 10, 60, 136, 23, 114}, {12, 80, 231, 208, 169, 191, 87, 195, 125, 38, 181, 47, 217, 197, 85, 219, 221, 245, 8, 96, 186, 107, 206, 33, 145, 130, 86, 207, 45, 193, 101, 134, 102, 146, 150, 166, 251, 64, 39, 185, 127, 62, 21, 252, 100, 138, 54, 117, 70, 15, 68, 23, 228, 196, 89}, {24, 93, 107, 129, 132, 252, 200, 18, 173, 3, 40, 231, 189, 158, 145, 25, 69, 54, 234, 5, 120, 52, 218, 191, 174, 43, 207, 90, 35, 15, 136, 92, 115, 220, 239, 125, 76, 238, 101, 17, 133, 228, 149, 121, 44, 135, 212, 47, 175, 51, 146, 49, 162, 139, 116}, {48, 105, 127, 248, 77, 241, 224, 247, 64, 156, 95, 182, 236, 170, 150, 162, 11, 205, 212, 94, 134, 133, 213, 110, 239, 250, 45, 35, 30, 26, 218, 99, 130, 69, 108, 143, 40, 211, 206, 132, 229, 7, 144, 2, 96, 210, 254, 237, 154, 255, 221, 243, 128, 37, 190}, {96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59, 85, 150, 89, 44, 38, 193, 15, 26, 169, 145, 100, 36, 1, 96, 185, 223, 59}, {192, 222, 182, 151, 114, 110, 155, 27, 143, 160, 177, 237, 82, 75, 89, 88, 152, 70, 240, 103, 21, 123, 224, 251, 116, 212, 101, 136, 218, 145, 200, 144, 8, 78, 190, 217, 204, 183, 87, 172, 216, 12, 105, 225, 59, 170, 98, 242, 250, 180, 10, 211, 31, 168, 255}, {157, 95, 217, 133, 230, 130, 18, 2, 39, 190, 175, 23, 209, 25, 36, 4, 78, 97, 67, 46, 191, 50, 72, 8, 156, 194, 134, 92, 99, 100, 144, 16, 37, 153, 17, 184, 198, 200, 61, 32, 74, 47, 34, 109, 145, 141, 122, 64, 148, 94, 68, 218, 63, 7, 244}};

/**
 * @brief Computes 2 * PARAM_DELTA syndromes
 *
 * @param[out] syndromes Array of size 2 * PARAM_DELTA receiving the computed syndromes
 * @param[in] cdw Array of size PARAM_N1 storing the received vector
 */
static void compute_syndromes(uint16_t *syndromes, uint8_t *cdw) {
    for (size_t i = 0; i < 2 * PARAM_DELTA; ++i) {
        for (size_t j = 1; j < PARAM_N1; ++j) {
            syndromes[i] ^= PQCLEAN_HQC192_CLEAN_gf_mul(cdw[j], alpha_ij_pow[i][j - 1]);
        }
        syndromes[i] ^= cdw[0];
    }
}
#endif


#if 1
#else
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

        dd = PQCLEAN_HQC192_CLEAN_gf_mul(d, PQCLEAN_HQC192_CLEAN_gf_inverse(d_p));

        for (i = 1; (i <= mu + 1) && (i <= PARAM_DELTA); ++i) {
            sigma[i] ^= PQCLEAN_HQC192_CLEAN_gf_mul(dd, X_sigma_p[i]);
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
            d ^= PQCLEAN_HQC192_CLEAN_gf_mul(sigma[i], syndromes[mu + 1 - i]);
        }
    }

    return deg_sigma;
}
#endif

static const uint32_t VandermoreMatrix_T[280] = {
    0x01010101, 0x01010101, 0x01010101, 0x01010101, 0x00000001, 
    0xad478e01, 0x1b366cd8, 0xfae9cf83, 0x2c58b07d, 0x00000016, 
    0x36d84701, 0x587de983, 0xf3eb8b16, 0x24907af5, 0x00000009, 
    0xcf36ad01, 0xfb8b2c7d, 0x8a243df5, 0x59f2c356, 0x000000a6, 
    0x7d83d801, 0x90f5eb16, 0xb2ef5609, 0x640ee0a6, 0x00000041, 
    0x2ce96c01, 0xac24f4eb, 0x1ca759ef, 0x9637ae64, 0x000000ff, 
    0x8b7d3601, 0xf25624f5, 0x576438a6, 0x91b3dbc4, 0x00000073, 
    0xfb581b01, 0x53f2ac90, 0xab37820e, 0x55d5c6b3, 0x00000054, 
    0xf5168301, 0x0ea6ef09, 0x7effc441, 0xa954aa73, 0x000000cc, 
    0x3d8bcf01, 0x82385956, 0xe491f1c4, 0x3b172129, 0x000000ce, 
    0x24ebe901, 0x3764a7ef, 0x4d7291ff, 0x1a7c33a9, 0x000000e2, 
    0x8af3fa01, 0xab571cb2, 0xda4de47e, 0xdf2281c5, 0x000000f0, 
    0x56f57d01, 0xb3c464a6, 0xc5a92973, 0x0f7f86ce, 0x000000a1, 
    0xc37ab001, 0xc6dbaee0, 0x813321aa, 0xb9cab186, 0x00000023, 
    0xf2905801, 0xd5b3370e, 0x227c1754, 0xc1d2ca7f, 0x000000c0, 
    0x59242c01, 0x55919664, 0xdf1a3ba9, 0x60c1b90f, 0x00000026, 
    0xa6091601, 0x5473ff41, 0xf0e2cecc, 0x26c023a1, 0x0000008e, 
    0xdd450b01, 0x4f92d7dc, 0x99d64493, 0x01984e0a, 0x0000000b, 
    0x38568b01, 0x172991c4, 0xba0fb6ce, 0x2c087535, 0x0000008a, 
    0x079bcb01, 0xec15e64b, 0x9fbe6b88, 0x24fa1d0c, 0x00000053, 
    0x64efeb01, 0x7ca972ff, 0x9ca00fe2, 0x59f46c26, 0x000000ae, 
    0x82f2fb01, 0xd01755b3, 0x8fc1617f, 0x64c38b08, 0x000000f1, 
    0x57b2f301, 0x22c54d7e, 0x879cbaf0, 0x967012cf, 0x000000b7, 
    0x6ea2f701, 0x71ed8463, 0x040346bc, 0x91a5b2fb, 0x0000002a, 
    0xc4a6f501, 0x7fcea973, 0xcf2635a1, 0x55f1078a, 0x00000066, 
    0x96a7f401, 0xfd1a2e72, 0xeb2060a0, 0xa9e63759, 0x00000067, 
    0xdbe07a01, 0xca8633aa, 0x126c7523, 0x3b9a7b07, 0x00000071, 
    0xf1383d01, 0x61b63b29, 0xf22ccd35, 0x1ab8736e, 0x00000078, 
    0xb30e9001, 0xd27f7c54, 0x70f408c0, 0xdfc79af1, 0x000000de, 
    0xfc8d4801, 0x14e76742, 0x41ac36c9, 0x0f0d6dbf, 0x0000009f, 
    0x91642401, 0xc10f1aa9, 0x96592c26, 0xb9df3b55, 0x00000060, 
    0xbf191201, 0x942f115c, 0xe51cf580, 0xc11e6821, 0x00000013, 
    0x73410901, 0xc0a1e2cc, 0xb7ae8a8e, 0x60de7166, 0x00000047, 
    0xe4578a01, 0x8fbadfc5, 0x2996f2cf, 0x2646e73e, 0x0000008b, 
    0x92dc4501, 0x980ad693, 0x4fd7dd0b, 0x014e9944, 0x00000045, 
    0x5537ac01, 0x74c1fd7c, 0x33e664f4, 0x2cb4a0df, 0x000000a7, 
    0x29c45601, 0x08350fce, 0x3e556e8a, 0x2440b578, 0x00000057, 
    0xa8312b01, 0xd8275e68, 0x8884dbf9, 0x59833061, 0x000000f6, 
    0x154b9b01, 0xfa0cbe88, 0x5b2efc53, 0x64f34c50, 0x000000d5, 
    0x21dbc301, 0x8b75b986, 0xe73b7307, 0x965608b5, 0x00000015, 
    0xa9ffef01, 0xf426a0e2, 0x5e6755ae, 0x91a7e960, 0x00000033, 
    0xb8f6f901, 0x093a05a3, 0xd2111562, 0x5541f72d, 0x000000bd, 
    0x17b3f201, 0xc308c17f, 0x46dfb8f1, 0xa9db5640, 0x000000b6, 
    0x66e57901, 0xa2ad6abb, 0x4afdc53f, 0x3b635336, 0x0000003c, 
    0xc57eb201, 0x70cf9cf0, 0x065e3eb7, 0x1aa4198b, 0x0000006f, 
    0x3b915901, 0x642c600f, 0x26b91a55, 0xdfa99624, 0x000000c1, 
    0xed63a201, 0xa5fb03bc, 0x1005d92a, 0x0f767ef2, 0x00000030, 
    0x3ed15101, 0x313db4c2, 0x1b6a7fda, 0xb9683938, 0x00000087, 
    0xce73a601, 0xf18a26a1, 0x8b607866, 0xc1b61557, 0x000000ad, 
    0xd0d55301, 0x7ec374d2, 0x48b42fc7, 0x60f085db, 0x000000cb, 
    0x1a72a701, 0xe65920a0, 0xef74b967, 0x26be7c91, 0x000000ac, 
    0x4492dd01, 0x92dd010a, 0xdd010a44, 0x010a4492, 0x000000dd, 
    0x86aae001, 0x9a076c23, 0x19e9b571, 0x2c4ae115, 0x000000a5, 
    0xd9a47001, 0x4282e977, 0x62eb27b1, 0x24c91e17, 0x0000007b, 
    0xb6293801, 0xb86e2c35, 0xb3248f78, 0x593aa1ed, 0x000000e4, 
    0xdf4d1c01, 0x3396eb9c, 0xe6ef265e, 0x646c051a, 0x00000084
};

// skipping the first column which are all equal to one
static const uint32_t VandermoreMatrix_T_56x16[224] = {
    0x1010101, 0x1010101, 0x1010101, 0x1010101, 
    0xd8ad478e, 0x831b366c, 0x7dfae9cf, 0x162c58b0, 
    0x8336d847, 0x16587de9, 0xf5f3eb8b, 0x924907a, 
    0x7dcf36ad, 0xf5fb8b2c, 0x568a243d, 0xa659f2c3, 
    0x167d83d8, 0x990f5eb, 0xa6b2ef56, 0x41640ee0, 
    0xeb2ce96c, 0xefac24f4, 0x641ca759, 0xff9637ae, 
    0xf58b7d36, 0xa6f25624, 0xc4576438, 0x7391b3db, 
    0x90fb581b, 0xe53f2ac, 0xb3ab3782, 0x5455d5c6, 
    0x9f51683, 0x410ea6ef, 0x737effc4, 0xcca954aa, 
    0x563d8bcf, 0xc4823859, 0x29e491f1, 0xce3b1721, 
    0xef24ebe9, 0xff3764a7, 0xa94d7291, 0xe21a7c33, 
    0xb28af3fa, 0x7eab571c, 0xc5da4de4, 0xf0df2281, 
    0xa656f57d, 0x73b3c464, 0xcec5a929, 0xa10f7f86, 
    0xe0c37ab0, 0xaac6dbae, 0x86813321, 0x23b9cab1, 
    0xef29058, 0x54d5b337, 0x7f227c17, 0xc0c1d2ca, 
    0x6459242c, 0xa9559196, 0xfdf1a3b, 0x2660c1b9, 
    0x41a60916, 0xcc5473ff, 0xa1f0e2ce, 0x8e26c023, 
    0xdcdd450b, 0x934f92d7, 0xa99d644, 0xb01984e, 
    0xc438568b, 0xce172991, 0x35ba0fb6, 0x8a2c0875, 
    0x4b079bcb, 0x88ec15e6, 0xc9fbe6b, 0x5324fa1d, 
    0xff64efeb, 0xe27ca972, 0x269ca00f, 0xae59f46c, 
    0xb382f2fb, 0x7fd01755, 0x88fc161, 0xf164c38b, 
    0x7e57b2f3, 0xf022c54d, 0xcf879cba, 0xb7967012, 
    0x636ea2f7, 0xbc71ed84, 0xfb040346, 0x2a91a5b2, 
    0x73c4a6f5, 0xa17fcea9, 0x8acf2635, 0x6655f107, 
    0x7296a7f4, 0xa0fd1a2e, 0x59eb2060, 0x67a9e637, 
    0xaadbe07a, 0x23ca8633, 0x7126c75, 0x713b9a7b, 
    0x29f1383d, 0x3561b63b, 0x6ef22ccd, 0x781ab873, 
    0x54b30e90, 0xc0d27f7c, 0xf170f408, 0xdedfc79a, 
    0x42fc8d48, 0xc914e767, 0xbf41ac36, 0x9f0f0d6d, 
    0xa9916424, 0x26c10f1a, 0x5596592c, 0x60b9df3b, 
    0x5cbf1912, 0x80942f11, 0x21e51cf5, 0x13c11e68, 
    0xcc734109, 0x8ec0a1e2, 0x66b7ae8a, 0x4760de71, 
    0xc5e4578a, 0xcf8fbadf, 0x3e2996f2, 0x8b2646e7, 
    0x9392dc45, 0xb980ad6, 0x444fd7dd, 0x45014e99, 
    0x7c5537ac, 0xf474c1fd, 0xdf33e664, 0xa72cb4a0, 
    0xce29c456, 0x8a08350f, 0x783e556e, 0x572440b5, 
    0x68a8312b, 0xf9d8275e, 0x618884db, 0xf6598330, 
    0x88154b9b, 0x53fa0cbe, 0x505b2efc, 0xd564f34c, 
    0x8621dbc3, 0x78b75b9, 0xb5e73b73, 0x15965608, 
    0xe2a9ffef, 0xaef426a0, 0x605e6755, 0x3391a7e9, 
    0xa3b8f6f9, 0x62093a05, 0x2dd21115, 0xbd5541f7, 
    0x7f17b3f2, 0xf1c308c1, 0x4046dfb8, 0xb6a9db56, 
    0xbb66e579, 0x3fa2ad6a, 0x364afdc5, 0x3c3b6353, 
    0xf0c57eb2, 0xb770cf9c, 0x8b065e3e, 0x6f1aa419, 
    0xf3b9159, 0x55642c60, 0x2426b91a, 0xc1dfa996, 
    0xbced63a2, 0x2aa5fb03, 0xf21005d9, 0x300f767e, 
    0xc23ed151, 0xda313db4, 0x381b6a7f, 0x87b96839, 
    0xa1ce73a6, 0x66f18a26, 0x578b6078, 0xadc1b615, 
    0xd2d0d553, 0xc77ec374, 0xdb48b42f, 0xcb60f085, 
    0xa01a72a7, 0x67e65920, 0x91ef74b9, 0xac26be7c, 
    0xa4492dd, 0x4492dd01, 0x92dd010a, 0xdd010a44, 
    0x2386aae0, 0x719a076c, 0x1519e9b5, 0xa52c4ae1, 
    0x77d9a470, 0xb14282e9, 0x1762eb27, 0x7b24c91e, 
    0x35b62938, 0x78b86e2c, 0xedb3248f, 0xe4593aa1, 
    0x9cdf4d1c, 0x5e3396eb, 0x1ae6ef26, 0x84646c05,
};

static const uint32_t VandermoreMatrix_T_56x16_rdx16[224] = {
    0xc280ce6b, 0x268a2cc5, 0xc862cead, 0x6c8cc82b, 0x26e0ec03, 0x8266a6cb, 0x8a8ee0cd, 0xc88a0a6d, 0x6488ec6f, 0x66e06827, 0xe666c0e3, 0x2e6e4e87, 0xca6620a9, 0xac8e4081, 0x0868c609, 0x80842cc3, 
    0x156bd6a8, 0x5e6ae217, 0xa9d9f5ab, 0xe66b4356, 0xd5f49e11, 0x9bdb07a1, 0x7a6b44a8, 0x6460dee5, 0x85276194, 0x53f81fdd, 0x9297c6ba, 0xb6d80798, 0xfb526ce2, 0xa86e1d7e, 0x691e0d4a, 0x4f6e4c67, 
    0xfea8e4d7, 0xeaa961fe, 0x9175c785, 0xa6ac36e3, 0x353cdb37, 0x1a7c4b98, 0xc0a293fa, 0x6faaaea3, 0xfdc60906, 0x5b37993d, 0xcf59b792, 0xa971681f, 0x9072103e, 0x03a727c4, 0x985ef8ed, 0xffa85560, 
    0x43d7b901, 0x35d74941, 0x6f2b817c, 0x5cd02b30, 0xc17a4fe2, 0xf52e416e, 0x7e1b0e5b, 0xc3d3c259, 0x00da92af, 0x127eb2c5, 0xd92889ad, 0x5322bffa, 0xf2cca13c, 0xef135179, 0x6ac45117, 0x32d008c5, 
    0xa801e1fe, 0x860a96ad, 0x54366d9c, 0x6f04c689, 0xfd971dd7, 0x413d285b, 0x8d83cf66, 0xf304b364, 0x3da297b8, 0xd59caef2, 0xf13fe8ac, 0x1032ca43, 0x063ace5b, 0xd98ce080, 0x00b02b14, 0x380871fc, 
    0xfffe58ff, 0xf5f611f1, 0x4b928be2, 0x5bf733fd, 0x87152e16, 0xba949e4f, 0x4dbd1128, 0xb7fecc54, 0x2c6239b5, 0x7412228c, 0x70015f66, 0xaa98e3ba, 0xc19a51e4, 0xd7b3094c, 0xa54b1629, 0x74f693b6, 
    0x57ffad29, 0x72f88554, 0xfb8e256a, 0x25f5a079, 0x10fad4eb, 0xb18c5bf0, 0x6fd5753d, 0x50f8082b, 0x8e738b0b, 0x02f0ab10, 0xddca51f5, 0x1c84dab6, 0x0aee25ff, 0xf5dba461, 0x416be1a7, 0x0ffa8556, 
    0x14294280, 0x4228d71c, 0xf4dfcc9b, 0x2d29a942, 0x055b356c, 0x45d1bfff, 0xb30d0c12, 0xd5204224, 0x88ec3e48, 0x535bed02, 0xad930a2c, 0x53d7464f, 0x1ce3f436, 0x330aefb7, 0xb923755f, 0x532d55dc, 
    0xc280da43, 0x278e2ecd, 0x6d5b0c39, 0x7d8a9d26, 0xff9ece7c, 0xd8546069, 0x798a6c44, 0xdc887d76, 0x7ec63a52, 0xfa961ff2, 0x12d963bf, 0x8b5ccfd8, 0xe5d32447, 0x97893f70, 0x01d08c42, 0xc3804bdf, 
    0xa94361fe, 0x9b4ca1ae, 0x40d2520b, 0xb5407990, 0xa4aef779, 0x08d92f4d, 0xe626f0a8, 0x564658be, 0x6fbe73dc, 0x41a9a7a9, 0xca9379ba, 0x87db180c, 0x1255a301, 0x6b2a08e6, 0x863c4b84, 0x674f0f5b, 
    0x7ffed13d, 0xfdf71d7b, 0x58d65aa3, 0xd3faa2fe, 0x3a908486, 0x80d4e958, 0xb39141a6, 0x36fa99da, 0xa9272d0d, 0xab9eab3f, 0x4f0593f0, 0x03d84c84, 0x5f5e16ed, 0x3e9221bf, 0x34ada48f, 0x60ffce35, 
    0x3c3d403c, 0xcb3d113e, 0x33ecebc7, 0xb93366c4, 0x94c91922, 0x33eafd3b, 0x25f9060a, 0x9c3da5b5, 0x2a9ff9da, 0x49c5d396, 0x004a85d8, 0x32ea9431, 0x42afb669, 0x5af2c324, 0xc48df67e, 0xc532199d, 
    0xbc3c8b7e, 0xc73f04b5, 0x5d063930, 0x743480c1, 0x04740518, 0xde0ab85e, 0xb77d46e1, 0x4f300f76, 0x2173f4c6, 0x43760608, 0xbe28fa3b, 0xe50f41d7, 0x4a232f05, 0x7477a3b0, 0x15600229, 0xfc38fb4d, 
    0x157ea800, 0x5870bc1a, 0xfc2ecbbc, 0x847e8e57, 0xcdf0f38c, 0xc2282bf1, 0x34ab30d8, 0x4b7fa885, 0x30d12fe7, 0xddf0dacd, 0xace318e3, 0x222da7c7, 0x31344c4b, 0x4da2f834, 0x0ad23d02, 0xb6744b45,
};


#if defined(_M4_ASM_)
void radix_16_matrix_vector_56x16(uint32_t *in, const uint32_t *mat, uint32_t *vec);
#else
void radix_16_matrix_vector_56x16(uint32_t *in, const uint32_t *mat, uint32_t *vec) {
    for(int i=0;i<(14);i++) { in[i] = radix_16_matrix_vector(&mat[16*i], vec, 16); }
}
#endif

static void compute_roots_new(uint8_t *error, uint8_t *sigma) {
    
    /*
    uint8_t * vande_u8 = (uint8_t *)VandermoreMatrix_T_56x16;
    for(int i = 0; i < PARAM_N1; i++) {
        uint8_t value = sigma[0];
        for(int j = 0; j < PARAM_DELTA; j++)
        {
            value ^= PQCLEAN_HQC192_CLEAN_gf_mul(sigma[j+1], vande_u8[i * (16) + j]);// 20 is because of zero padding in VandermoreMatrix_T
        }
        error[i] = (value == 0);
    }
    */

    uint32_t _sigma_rdx16[4];
    uint8_t sigma_new[16] = {0};
    for(int i=0; i<16; i++) { sigma_new[i] = sigma[i+1]; }
    memcpy( _sigma_rdx16 , sigma_new , 16 );
    for(int i=0;i<4;i++) { _sigma_rdx16[i] = rdx16_from_bitseq( _sigma_rdx16[i] ); }
    uint32_t sigma_rdx16[16];
    for(int i=0; i<4; i++) {
        sigma_rdx16[4*i] = _sigma_rdx16[i]&_NORM_RDX16_;
        sigma_rdx16[4*i+1] = (_sigma_rdx16[i]>>2)&_NORM_RDX16_;
        sigma_rdx16[4*i+2] = (_sigma_rdx16[i]>>1)&_NORM_RDX16_;
        sigma_rdx16[4*i+3] = (_sigma_rdx16[i]>>3)&_NORM_RDX16_;
    }
    uint32_t sigma_0_rdx16 = rdx16_from_bitseq((uint32_t)sigma[0]);
    uint32_t sigma_0_rdx16_2 = sigma_0_rdx16<<2;
    uint32_t sigma_0_rdx16_1 = sigma_0_rdx16<<1;
    uint32_t sigma_0_rdx16_3 = sigma_0_rdx16<<3;
        
    uint32_t product[14] = {0};
    radix_16_matrix_vector_56x16(&product[0], &VandermoreMatrix_T_56x16_rdx16[0], sigma_rdx16);
    for(int i=0; i<14; i++) { 
        error[4*i] = (((product[i]&0x11111111)^(sigma_0_rdx16)) == 0);
        error[4*i+1] = (((product[i]&0x44444444)^(sigma_0_rdx16_2)) == 0);
        error[4*i+2] = (((product[i]&0x22222222)^(sigma_0_rdx16_1)) == 0);
        error[4*i+3] = (((product[i]&0x88888888)^(sigma_0_rdx16_3)) == 0);
    }
}

/**
 * @brief Computes the error polynomial error from the error locator polynomial sigma
 *
 * See function PQCLEAN_HQC192_CLEAN_fft for more details.
 *
 * @param[out] error Array of 2^PARAM_M elements receiving the error polynomial
 * @param[out] error_compact Array of PARAM_DELTA + PARAM_N1 elements receiving a compact representation of the vector error
 * @param[in] sigma Array of 2^PARAM_FFT elements storing the error locator polynomial
 */
static void compute_roots(uint8_t *error, uint16_t *sigma) {
    uint16_t w[1 << PARAM_M] = {0};

    PQCLEAN_HQC192_CLEAN_fft(w, sigma, PARAM_DELTA + 1);
    PQCLEAN_HQC192_CLEAN_fft_retrieve_error_poly(error, w);
}


#if 1
#else
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
            z[i] ^= mask & PQCLEAN_HQC192_CLEAN_gf_mul(sigma[j], syndromes[i - j - 1]);
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
        inverse = PQCLEAN_HQC192_CLEAN_gf_inverse(beta_j[i]);
        inverse_power_j = 1;

        for (size_t j = 1; j <= PARAM_DELTA; ++j) {
            inverse_power_j = PQCLEAN_HQC192_CLEAN_gf_mul(inverse_power_j, inverse);
            tmp1 ^= PQCLEAN_HQC192_CLEAN_gf_mul(inverse_power_j, z[j]);
        }
        for (size_t k = 1; k < PARAM_DELTA; ++k) {
            tmp2 = PQCLEAN_HQC192_CLEAN_gf_mul(tmp2, (1 ^ PQCLEAN_HQC192_CLEAN_gf_mul(inverse, beta_j[(i + k) % PARAM_DELTA])));
        }
        mask1 = (uint16_t) (((int16_t) i - delta_real_value) >> 15); // i < delta_real_value
        e_j[i] = mask1 & PQCLEAN_HQC192_CLEAN_gf_mul(tmp1, PQCLEAN_HQC192_CLEAN_gf_inverse(tmp2));
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

#if 1
/**
 * @brief Correct the errors
 *
 * @param[out] cdw Array of PARAM_N1 elements receiving the corrected vector
 * @param[in] error Array of the error vector
 * @param[in] error_values Array of PARAM_DELTA elements storing the error values
 */
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
void PQCLEAN_HQC192_CLEAN_reed_solomon_decode(uint8_t *msg, uint8_t *cdw) {
#if 1
    uint8_t syndromes_new[2 * PARAM_DELTA] = {0};
    uint8_t sigma_new[1 << PARAM_FFT] = {0};
    // uint16_t sigma_new_16[1 << PARAM_FFT] = {0};
    uint8_t omega[1 << PARAM_FFT] = {0};
    uint8_t error_values_new[PARAM_N1] = {0};
    uint8_t error_new[1 << PARAM_M] = {0};
    // uint8_t error_new2[1 << PARAM_M] = {0};
    // uint8_t mask = 0;

    // Calculate the 2*PARAM_DELTA syndromes
    compute_syndromes(syndromes_new, cdw);
    // for (size_t i = 0; i < 2 * PARAM_DELTA; i++) {
    //     syndromes_new[i] = (uint8_t)syndromes[i];
    // }


    // Compute the error locator polynomial sigma
    // Sigma's degree is at most PARAM_DELTA but the FFT requires the extra room 
    compute_elp_divstep(omega, sigma_new, syndromes_new);
    // for(int i = 0; i < (1 << PARAM_FFT); i++)
    // {
    //     sigma_new_16[i] = sigma_new[i];
    //     mask ^= sigma_new[i];
    // }
    // Compute the error polynomial error
    // compute_roots(error_new, sigma_new_16);
    compute_roots_new(error_new, sigma_new);


    compute_error_values_new(error_values_new, omega, sigma_new, error_new);

    // Correct the errors
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
