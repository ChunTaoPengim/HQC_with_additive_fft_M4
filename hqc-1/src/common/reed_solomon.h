#ifndef REED_SOLOMON_H
#define REED_SOLOMON_H


/**
 * @file reed_solomon.h
 * @brief Header file of reed_solomon.c
 */

#include <stdint.h>

void PQCLEAN_HQC128_CLEAN_reed_solomon_encode(uint8_t *cdw, const uint8_t *msg);

void PQCLEAN_HQC128_CLEAN_reed_solomon_decode(uint8_t *msg, uint8_t *cdw);


#endif
