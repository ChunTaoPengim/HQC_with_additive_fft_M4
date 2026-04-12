#ifndef _BC_1_H_
#define _BC_1_H_

#include "stdint.h"


//
// libaray for basis conversion
// computation unit: 1 bit
//


void bc_1_256( void *poly , unsigned n_256bit );

void ibc_1_256( void *poly , unsigned n_256bit );


/////////////////////////////////////////

// n_byte >= 32
void bc_1( void * poly , unsigned n_byte );

void ibc_1( void * poly , unsigned n_byte );


/// for HQC ////

void bc_1_16384_plus_2048(uint32_t *poly);
void bc_1_32768_plus_4096(uint32_t *poly);

static inline void ibc_1_32768(uint32_t *poly) { ibc_1(poly,4096); }
static inline void ibc_1_65536(uint32_t *poly) { ibc_1(poly,8192); }
static inline void  bc_1_65536(uint32_t *poly) { bc_1(poly,8192); }
static inline void ibc_1_131072(uint32_t *poly) { ibc_1(poly,16384); }




#endif
