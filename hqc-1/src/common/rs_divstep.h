#ifndef RS_DIVSTEP_H
#define RS_DIVSTEP_H


/**
 * @file rs_divstep.h
 * @brief Header file of rs_divstep.c
 */

#include <stdint.h>

void compute_elp_divstep(uint8_t * omega, uint8_t *sigma, uint8_t *syndromes);
void compute_error_values_new(uint8_t *error_values, const uint8_t *omega, const uint8_t *sigma, const uint8_t *error);


#endif