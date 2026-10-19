#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "hal.h"
#include "parameters.h"
#include "reed_solomon.h"

#if PARAM_N1 == 46
#define RS_ENCODE PQCLEAN_HQC128_CLEAN_reed_solomon_encode
#define RS_DECODE PQCLEAN_HQC128_CLEAN_reed_solomon_decode
#elif PARAM_N1 == 56
#define RS_ENCODE PQCLEAN_HQC192_CLEAN_reed_solomon_encode
#define RS_DECODE PQCLEAN_HQC192_CLEAN_reed_solomon_decode
#else
#define RS_ENCODE PQCLEAN_HQC256_CLEAN_reed_solomon_encode
#define RS_DECODE PQCLEAN_HQC256_CLEAN_reed_solomon_decode
#endif

// Each trial encodes a random message, adds NUM_ERRORS symbol errors at
// distinct random positions, and checks that decoding returns the message.
#define NUM_ERRORS (5)
#define TEST_RUN (100)

static uint32_t rng_state = 0x12345678u;

static uint32_t rng_next(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static int decode_with_errors(int num_errors) {
    uint8_t msg[PARAM_K];
    uint8_t decoded[PARAM_K];
    uint8_t cdw[PARAM_N1] = {0};
    uint8_t flipped[PARAM_N1] = {0};

    for (int i = 0; i < PARAM_K; i++) {
        msg[i] = (uint8_t)rng_next();
    }
    RS_ENCODE(cdw, msg);

    for (int e = 0; e < num_errors; e++) {
        int pos;
        uint8_t err;
        do {
            pos = rng_next() % PARAM_N1;
        } while (flipped[pos]);
        flipped[pos] = 1;
        do {
            err = (uint8_t)rng_next();
        } while (err == 0);
        cdw[pos] ^= err;
    }

    RS_DECODE(decoded, cdw);
    return memcmp(msg, decoded, PARAM_K) == 0;
}

int main(void) {
    char outstr[128];
    int passed = 0;

    hal_setup(CLOCK_BENCHMARK);

    for (int i = 0; i < TEST_RUN; i++) {
        passed += decode_with_errors(NUM_ERRORS);
    }

    sprintf(outstr, "[%d/%d] RS decoding with %d symbol errors passed\n", passed, TEST_RUN, NUM_ERRORS);
    hal_send_str(outstr);
    if (passed == TEST_RUN) {
        hal_send_str("\n*** ALL GOOD ***\n");
    } else {
        hal_send_str("\n*** TEST FAILED ***\n");
    }

    return 0;
}
