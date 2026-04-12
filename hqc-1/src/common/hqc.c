#include "code.h"
#include "gf2x.h"
#include "hqc.h"
#include "parameters.h"
#include "parsing.h"
#include "symmetric.h"
#include "crypto_memset.h"
#include "vector.h"
#include "gfmul_fft.h"
#include <stdint.h>
#include <string.h>
/**
 * @file hqc.c
 * @brief Implementation of hqc.h
 */



/**
 * @brief Keygen of the HQC_PKE IND_CPA scheme
 *
 * The public key is composed of the syndrome <b>s</b> as well as the <b>seed</b> used to generate the vector <b>h</b>.
 *
 * The secret key is composed of the <b>seed</b> used to generate vectors <b>x</b> and  <b>y</b>.
 * As a technicality, the public key is appended to the secret key in order to respect NIST API.
 *
 * @param[out] pk String containing the public key
 * @param[out] sk String containing the secret key
 */
void PQCLEAN_HQC128_CLEAN_hqc_pke_keygen(uint8_t *ek_pke, uint8_t *dk_pke, uint8_t *seed) {
    uint8_t keypair_seed[2 * SEED_BYTES] = {0};
    uint8_t *seed_dk = keypair_seed;
    uint8_t *seed_ek = keypair_seed + SEED_BYTES;
    shake256_xof_ctx dk_xof_ctx = {0};
    shake256_xof_ctx ek_xof_ctx = {0};

    uint64_t x[VEC_N_SIZE_64] = {0};
    uint64_t y[VEC_N_SIZE_64] = {0};
    uint64_t h[VEC_N_SIZE_64] = {0};
    uint64_t s[VEC_N_SIZE_64] = {0};

    // Derive keypair seeds
    hash_i(keypair_seed, seed);

    // Compute decryption key
    xof_init(&dk_xof_ctx, seed_dk, SEED_BYTES);
    vect_sample_fixed_weight1(&dk_xof_ctx, y, PARAM_OMEGA);
    vect_sample_fixed_weight1(&dk_xof_ctx, x, PARAM_OMEGA);


    // Compute encryption key
    xof_init(&ek_xof_ctx, seed_ek, SEED_BYTES);
    vect_set_random(&ek_xof_ctx, h);
    PQCLEAN_HQC128_CLEAN_vect_mul(s, y, h);
    vect_add(s, x, s, VEC_N_SIZE_64);

    // Parse encryption key to string
    memcpy(ek_pke, seed_ek, SEED_BYTES);
    memcpy(ek_pke + SEED_BYTES, s, VEC_N_SIZE_BYTES);

    // Parse decryption key to string
    memcpy(dk_pke, seed_dk, SEED_BYTES);

    // Zeroize sensitive data
    memset_zero(keypair_seed, sizeof keypair_seed);
    memset_zero(x, sizeof x);
    memset_zero(y, sizeof y);
    memset_zero(&dk_xof_ctx, sizeof dk_xof_ctx);


}



/**
 * @brief Encryption of the HQC_PKE IND_CPA scheme
 *
 * The cihertext is composed of vectors <b>u</b> and <b>v</b>.
 *
 * @param[out] u Vector u (first part of the ciphertext)
 * @param[out] v Vector v (second part of the ciphertext)
 * @param[in] m Vector representing the message to encrypt
 * @param[in] theta Seed used to derive randomness required for encryption
 * @param[in] pk String containing the public key
 */
void PQCLEAN_HQC128_CLEAN_hqc_pke_encrypt(ciphertext_pke_t *c_pke, const uint8_t *ek_pke, const uint64_t *m, const uint8_t *theta) {
   shake256_xof_ctx theta_xof_ctx = {0};
    uint64_t h[VEC_N_SIZE_64] = {0};
    uint64_t s[VEC_N_SIZE_64] = {0};
    uint64_t r1[VEC_N_SIZE_64] = {0};
    uint64_t r2[VEC_N_SIZE_64] = {0};
    uint64_t e[VEC_N_SIZE_64] = {0};
    uint64_t tmp[VEC_N_SIZE_64] = {0};


    // Initialize Xof using theta
    xof_init(&theta_xof_ctx, theta, SEED_BYTES);

    // Retrieve h and s from public key
    hqc_ek_pke_from_string(h, s, ek_pke);

    // Generate re, e and r1
    vect_sample_fixed_weight2(&theta_xof_ctx, r2, PARAM_OMEGA_R);
    vect_sample_fixed_weight2(&theta_xof_ctx, e, PARAM_OMEGA_E);
    vect_sample_fixed_weight2(&theta_xof_ctx, r1, PARAM_OMEGA_R);


    // Compute u = r1 + r2.h
    // Compute v = m.G + s.r2 + e
    PQCLEAN_HQC128_CLEAN_vect_mul_2(c_pke->u, tmp, r2, h, s);
    vect_add(c_pke->u, r1, c_pke->u, VEC_N_SIZE_64);

    // Compute v = C.encode(m)
    PQCLEAN_HQC128_CLEAN_code_encode(c_pke->v, m);

    // Compute v = C.encode(m) + Truncate(s.r2 + e)
    vect_add(tmp, e, tmp, VEC_N_SIZE_64);
    vect_truncate(tmp);
    vect_add(c_pke->v, c_pke->v, tmp, VEC_N1N2_SIZE_64);

    // Zeroize sensitive data
    memset_zero(r1, sizeof r1);
    memset_zero(r2, sizeof r2);
    memset_zero(e, sizeof e);
    memset_zero(tmp, sizeof tmp);
    memset_zero(&theta_xof_ctx, sizeof theta_xof_ctx);


}



/**
 * @brief Decryption of the HQC_PKE IND_CPA scheme
 *
 * @param[out] m Vector representing the decrypted message
 * @param[in] u Vector u (first part of the ciphertext)
 * @param[in] v Vector v (second part of the ciphertext)
 * @param[in] sk String containing the secret key
 * @returns 0
 */
uint8_t PQCLEAN_HQC128_CLEAN_hqc_pke_decrypt(uint64_t *m, const uint8_t *dk_pke, const ciphertext_pke_t *c_pke) {
    uint64_t y[VEC_N_SIZE_64] = {0};
    uint64_t tmp1[VEC_N_SIZE_64] = {0};
    uint64_t tmp2[VEC_N_SIZE_64] = {0};

    // Parse decryption key dk_pke
    hqc_dk_pke_from_string(y, dk_pke);

    // Compute u.y
    PQCLEAN_HQC128_CLEAN_vect_mul(tmp1, y, c_pke->u);
    // Truncate(u.y)
    vect_truncate(tmp1);
    // Compute v - Truncate(u.y)
    vect_add(tmp2, c_pke->v, tmp1, VEC_N1N2_SIZE_64);


    // Compute plaintext m
    PQCLEAN_HQC128_CLEAN_code_decode(m, tmp2);

    // Zeroize sensitive data
    memset_zero(y, sizeof y);
    memset_zero(tmp1, sizeof tmp1);
    memset_zero(tmp2, sizeof tmp2);

    return 0;
}
