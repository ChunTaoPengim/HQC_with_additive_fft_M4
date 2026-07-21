// SPDX-License-Identifier: Apache-2.0 or CC0-1.0
#include "api.h"
#include "hal.h"
#include "sendfn.h"

#include <stdint.h>
#include <string.h>

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_BYTES           NAMESPACE(CRYPTO_BYTES)
#define MUPQ_CRYPTO_PUBLICKEYBYTES  NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES  NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_CIPHERTEXTBYTES NAMESPACE(CRYPTO_CIPHERTEXTBYTES)
#define MUPQ_CRYPTO_ALGNAME NAMESPACE(CRYPTO_ALGNAME)

#define MUPQ_crypto_kem_keypair NAMESPACE(crypto_kem_keypair)
#define MUPQ_crypto_kem_enc NAMESPACE(crypto_kem_enc)
#define MUPQ_crypto_kem_dec NAMESPACE(crypto_kem_dec)
#define MUPQ_mul_size NAMESPACE(mul_size)
#define MUPQ_vect_mul NAMESPACE(PQCLEAN_HQC_CLEAN_vect_mul)
#define MUPQ_VEC_K_CODE NAMESPACE(VEC_K_CODE)
#define MUPQ_VEC_N1N2_code NAMESPACE(VEC_N1N2_code)
#define MUPQ_PQCLEAN_HQC_encode NAMESPACE(PQCLEAN_HQC_encode)
#define MUPQ_PQCLEAN_HQC_decode NAMESPACE(PQCLEAN_HQC_decode)

#define printcycles(S, U) send_unsignedll((S), (U))

int main(void)
{
  unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
  unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[MUPQ_CRYPTO_CIPHERTEXTBYTES];
  unsigned long long t0, t1;
  int i;

  hal_setup(CLOCK_BENCHMARK);

  hal_send_str("==========================");

  for(i=0;i<MUPQ_ITERATIONS; i++)
  {
    // Key-pair generation
    t0 = hal_get_time();
    MUPQ_crypto_kem_keypair(pk, sk);
    t1 = hal_get_time();
    printcycles("keypair cycles:", t1-t0);

    // Encapsulation
    t0 = hal_get_time();
    MUPQ_crypto_kem_enc(ct, key_a, pk);
    t1 = hal_get_time();
    printcycles("encaps cycles:", t1-t0);
    // Decapsulation
    t0 = hal_get_time();
    MUPQ_crypto_kem_dec(key_b, ct, sk);
    t1 = hal_get_time();
    printcycles("decaps cycles:", t1-t0);

    if (memcmp(key_a, key_b, MUPQ_CRYPTO_BYTES)) {
      hal_send_str("ERROR KEYS\n");
    }
    else {
      hal_send_str("OK KEYS\n");
    }
    hal_send_str("+");

    // Multiplication
    hal_send_str("test mul:");
    uint64_t a[MUPQ_mul_size];
    uint64_t b[MUPQ_mul_size];
    uint64_t o[MUPQ_mul_size];
    memset(a, 0xAA, sizeof(a));
    memset(b, 0xBB, sizeof(b));
    t0 = hal_get_time();
    MUPQ_vect_mul(o, a, b);
    t1 = hal_get_time();
    printcycles("mul cycles:", t1 - t0);

    // Encoding
    hal_send_str("test enc:");
    uint8_t msg[MUPQ_VEC_K_CODE];
    uint64_t codeword[MUPQ_VEC_N1N2_code];
    memset(msg, 0xCC, sizeof(msg));
    t0 = hal_get_time();
    MUPQ_PQCLEAN_HQC_encode(codeword, msg);
    t1 = hal_get_time();
    printcycles("enc cycles:", t1 - t0);

    // Decoding
    hal_send_str("test dec:");
    uint8_t rec_msg[MUPQ_VEC_K_CODE];
    t0 = hal_get_time();
    MUPQ_PQCLEAN_HQC_decode(rec_msg, codeword);
    t1 = hal_get_time();
    printcycles("dec cycles:", t1 - t0);
  }

  hal_send_str("#");

  return 0;
}
