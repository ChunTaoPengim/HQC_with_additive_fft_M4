
#include "stdio.h"

#include "parameters.h"

#include "gf.h"


#if defined(HQC128)
#define gf_mul         PQCLEAN_HQC128_CLEAN_gf_mul
#define gf_square      PQCLEAN_HQC128_CLEAN_gf_square
#define gf_inverse     PQCLEAN_HQC128_CLEAN_gf_inverse
#elif defined(HQC192)
#define gf_mul         PQCLEAN_HQC192_CLEAN_gf_mul
#define gf_square      PQCLEAN_HQC192_CLEAN_gf_square
#define gf_inverse     PQCLEAN_HQC192_CLEAN_gf_inverse
#elif defined(HQC256)
#define gf_mul         PQCLEAN_HQC256_CLEAN_gf_mul
#define gf_square      PQCLEAN_HQC256_CLEAN_gf_square
#define gf_inverse     PQCLEAN_HQC256_CLEAN_gf_inverse
#else
    #error "Unsupported SECURITYLEVEL"
#endif


int main(void)
{

    for (int a = 0; a < 256; a++) {
        for (int b = 0; b < 256; b++) {
            uint16_t res1 = gf_mul((uint16_t)a, (uint16_t)b);
            uint16_t res2 = 0;
            // naive carryless multiplication
            for (int i = 0; i < 8; i++) {
                if ((b >> i) & 1) {
                    res2 ^= ((uint16_t)a) << i;
                }
            }
            // reduction modulo x^8 + x^4 + x^3 + x^2 + 1
            for (int i = 15; i >= 8; i--) {
                if ((res2 >> i) & 1) {
                    res2 ^= 0x11D << (i - 8);
                }
            }
            if (res1 != res2) {
                printf("Mismatch for a=%x, b=%x: res1=%x, res2=%x\n", a, b, res1, res2);
                return 1;
            }
        }
    }
    printf("All MUL tests passed!\n");


    for (int a = 0; a < 256; a++) {
        do {
            int b = a;
            uint16_t res1 = gf_square((uint16_t)a);
            uint16_t res2 = 0;
            // naive carryless multiplication
            for (int i = 0; i < 8; i++) {
                if ((b >> i) & 1) {
                    res2 ^= ((uint16_t)a) << i;
                }
            }
            // reduction modulo x^8 + x^4 + x^3 + x^2 + 1
            for (int i = 15; i >= 8; i--) {
                if ((res2 >> i) & 1) {
                    res2 ^= 0x11D << (i - 8);
                }
            }
            if (res1 != res2) {
                printf("Mismatch for a=%x, b=%x: res1=%x, res2=%x\n", a, b, res1, res2);
                return 1;
            }
        } while(0);
    }
    printf("All SQU tests passed!\n");


    for (int a = 0; a < 256; a++) {
        do {
            int b = gf_inverse((uint16_t)a);
            uint16_t res1 = gf_mul((uint16_t)a, (uint16_t)b);
            uint16_t res2 = 0;
            // naive carryless multiplication
            for (int i = 0; i < 8; i++) {
                if ((b >> i) & 1) {
                    res2 ^= ((uint16_t)a) << i;
                }
            }
            // reduction modulo x^8 + x^4 + x^3 + x^2 + 1
            for (int i = 15; i >= 8; i--) {
                if ((res2 >> i) & 1) {
                    res2 ^= 0x11D << (i - 8);
                }
            }
            if ( a == 0) {
                printf("Testing inverse of zero, got %x\n", b);
            } else if (res1 != 1) {
                printf("Mismatch for a=%x, b=%x: res1=%x, res2=1\n", a, b, res1);
                return 1;
            }
            if (res1 != res2) {
                printf("Mismatch for a=%x, b=%x: res1=%x, res2=%x\n", a, b, res1, res2);
                return 1;
            }
        } while(0);
    }
    printf("All INV tests passed!\n");


    return 0;

}

