#ifndef _BS_RS_ENCODING_H_
#define _BS_RS_ENCODING_H_

#include "stdint.h"


void rs_encode_46_16( uint8_t * codeword , const uint8_t * mesg );

void rs_encode_56_24( uint8_t * codeword , const uint8_t * mesg );

void rs_encode_90_32( uint8_t * codeword , const uint8_t * mesg );


#endif

