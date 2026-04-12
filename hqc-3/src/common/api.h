#ifndef PQCLEAN_HQC192_CLEAN_API_H
#define PQCLEAN_HQC192_CLEAN_API_H
/**
 * @file api.h
 * @brief NIST KEM API used by the HQC_KEM IND-CCA2 scheme
 */

#include <stdint.h>

#define PQCLEAN_HQC192_CLEAN_CRYPTO_ALGNAME                      "HQC-192"

#define CRYPTO_SECRETKEYBYTES               4602
#define CRYPTO_PUBLICKEYBYTES               4514
#define CRYPTO_BYTES                        32
#define CRYPTO_CIPHERTEXTBYTES              8978

// As a technicality, the public key is appended to the secret key in order to respect the NIST API.
// Without this constraint, PQCLEAN_HQC192_CLEAN_CRYPTO_SECRETKEYBYTES would be defined as 32

int crypto_kem_keypair(uint8_t *pk, uint8_t *sk);

int crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk);

int crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);


#endif
