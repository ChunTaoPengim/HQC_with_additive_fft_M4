#ifndef _RS_SYNDROME_H_
#define _RS_SYNDROME_H_

#include "stdint.h"

#define SYNDROME_T uint8_t

void rs_syndrome_n46_r30( SYNDROME_T * syndromes , const uint8_t * cdw );

void rs_syndrome_n56_r32( SYNDROME_T * syndromes , const uint8_t * cdw );

void rs_syndrome_n90_r58( SYNDROME_T * syndromes , const uint8_t * cdw );


#endif

