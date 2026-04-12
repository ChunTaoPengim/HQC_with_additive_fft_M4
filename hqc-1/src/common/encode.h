// Implemented by Ming-Shing Chen, Tung Chou and Markus Krausz.
// public domain

#ifndef _ENCODE_H_
#define _ENCODE_H_

#include "stdint.h"

#include "gfv_tower.h"  // for sto_t


void encode_to_gft_full_length(uint32_t * out , const uint32_t * in );
void decode_from_gft_full_length(uint32_t * out , const uint32_t * in );


#endif // _ENCODE_H_
