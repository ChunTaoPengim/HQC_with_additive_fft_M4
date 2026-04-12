#ifndef HQC_H
#define HQC_H


/**
 * @file hqc.h
 * @brief Functions of the HQC_PKE IND_CPA scheme
 */

#include <stdint.h>
#include "data_structures.h"

void PQCLEAN_HQC128_CLEAN_hqc_pke_keygen(uint8_t *ek_pke, uint8_t *dk_pke, uint8_t *seed);

void PQCLEAN_HQC128_CLEAN_hqc_pke_encrypt(ciphertext_pke_t *c_pke, const uint8_t *ek_pke, const uint64_t *m, const uint8_t *theta);

uint8_t PQCLEAN_HQC128_CLEAN_hqc_pke_decrypt(uint64_t *m, const uint8_t *dk_pke, const ciphertext_pke_t *c_pke);


#endif
