#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "hal.h"
#include "api.h"
#include "parameters.h"
#include "gfmul_fft.h"
#include "radix_mul.h"


static int vec_equal(unsigned char *a, unsigned char *b, int len) {
    for(int i = 0 ; i < len ; ++i) {
        if(a[i] != b[i]) return 0;
    }
    return 1;
}


#define TEST_RUN (1000)


int main(void)
{
    
	printf("N: %d   ", PARAM_N);
	printf("N1: %d   ", PARAM_N1);
	printf("N2: %d   ", PARAM_N2);
	printf("OMEGA: %d   ", PARAM_OMEGA);
	printf("OMEGA_R: %d   ", PARAM_OMEGA_R);
    printf("\n");

    unsigned char pk[PUBLIC_KEY_BYTES];
	unsigned char sk[SECRET_KEY_BYTES];
	unsigned char ct[CIPHERTEXT_BYTES];
	unsigned char key1[SHARED_SECRET_BYTES];
	unsigned char key2[SHARED_SECRET_BYTES];
    hal_setup(CLOCK_BENCHMARK);

    int passed = 1;
    for(int i=0;i<TEST_RUN;i++) {
        
    	crypto_kem_keypair(pk, sk);
        crypto_kem_enc(ct, key1, pk);

	    crypto_kem_dec(key2, ct, sk);

        if(!vec_equal(key1, key2, SHARED_SECRET_BYTES)) {
            printf("[%d] Error: key1 != key2\n", i );
            passed = 0;
            hal_send_str("\n*** TEST FAILED ***\n");
            return -1;
        }


    }

    uint64_t cycles;
    char cycles_str[64];
    // benchmark keygen
    cycles = hal_get_time();
    crypto_kem_keypair(pk, sk);
    cycles = hal_get_time() - cycles;
#ifdef MPS2_AN386
    (void)cycles;
    sprintf(cycles_str, "[cycle counts not meaningful in qemu emulation]\n");
#else
    sprintf(cycles_str, "%llu\n", cycles);
#endif
    hal_send_str(cycles_str);

    // benchmark encap
    cycles = hal_get_time();
    crypto_kem_enc(ct, key1, pk);
    cycles = hal_get_time() - cycles;
#ifdef MPS2_AN386
    (void)cycles;
    sprintf(cycles_str, "[cycle counts not meaningful in qemu emulation]\n");
#else
    sprintf(cycles_str, "%llu\n", cycles);
#endif
    hal_send_str(cycles_str);

    // benchmark decap
    cycles = hal_get_time();
    crypto_kem_dec(key2, ct, sk);
    cycles = hal_get_time() - cycles;
#ifdef MPS2_AN386
    (void)cycles;
    sprintf(cycles_str, "[cycle counts not meaningful in qemu emulation]\n");
#else
    sprintf(cycles_str, "%llu\n", cycles);
#endif
    hal_send_str(cycles_str);



    if(passed) {
        char outstr[128];
        sprintf(outstr, "[%d/%d] KEM functional tests passed\n", TEST_RUN, TEST_RUN);
        hal_send_str(outstr);
        hal_send_str("\n*** ALL GOOD ***\n");
    }

    return 0;
}