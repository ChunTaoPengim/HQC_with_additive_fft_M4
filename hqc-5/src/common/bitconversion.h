// Implemented by Myeonghoon Lee and Jihoon Jang
// public domain

#ifndef _BITCONVERSION_H_
#define _BITCONVERSION_H_

#include <stdint.h>


void bc_1_16384(uint32_t *poly);
void ibc_1_16384(uint32_t *poly);
void bc_1_32768(uint32_t *poly);
void ibc_1_32768(uint32_t *poly);
void bc_1_65536(uint32_t *poly);
void ibc_1_65536(uint32_t *poly);
void bc_1_131072(uint32_t *poly);
void ibc_1_131072(uint32_t *poly);

#endif
