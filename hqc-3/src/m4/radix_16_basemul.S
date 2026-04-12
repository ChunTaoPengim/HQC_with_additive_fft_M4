// Implemented by Jihoon Jang and Myeonghoon Lee
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global radix_16_basemul_asm
.type radix_16_basemul_asm, %function
.align 2
radix_16_basemul_asm:

stmdb   sp!, { r4-r11, lr }
@ mask32 = 0x11111111
mov.w	r3, #286331153
@ a0 = (a)&mask32;
and.w	r4, r3, r1
@ b0 = (b)&mask32;
and.w	r6, r3, r2
@ b1 = (b >> 2) & mask32
and.w	r7, r3, r2, lsr #2
@ b2 = (b >> 1) & mask32
and.w	r8, r3, r2, lsr #1
@ b3 = (b >> 3) & mask32
and.w	r9, r3, r2, lsr #3
@ 0, [out[0], tmp[1]] = a0 * b0
umull	r2, r10, r4, r6
@ out[0] &= mask32
and.w	r2, r3, r2
@ 1, [tmp[2], tmp[3]] = a0 * b2
umull	r11, r12, r4, r8
@ 2, [tmp[1], tmp[2]] += a0 * b1
umlal	r10, r11, r4, r7
@ tmp[1] &= mask32
and.w	r10, r3, r10
@ a1 = (a >> 2) & mask32
and.w	r5, r3, r1, lsr #2
@ 3, [tmp[1], tmp[2]] += a1 * b0
umlal	r10, r11, r5, r6
@ tmp[1] &= mask32
and.w	r10, r3, r10
@ tmp[2] &= mask32
and.w	r11, r3, r11
@ out[0] ^= (tmp[1]<<2)
eor.w	r2, r2, r10, lsl #2
@ 4, [tmp[4], tmp[5]] = a1 * b3
umull	r14, r10, r5, r9
@ 5, [tmp[3], tmp[4]] += a0 * b3
umlal	r12, r14, r4, r9
@ 6, [tmp[2], tmp[3]] += a1 * b1
umlal	r11, r12, r5, r7
@ tmp[2] &= mask32
and.w	r11, r3, r11
@ tmp[3] &= mask32
and.w	r12, r3, r12
@ 7, [tmp[3], tmp[4]] += a1 * b2
umlal	r12, r14, r5, r8
@ tmp[4] &= mask32
and.w	r14, r3, r14
@ a2 = (a >> 1) & mask32
and.w	r5, r3, r1, lsr #1
@ 8, [tmp[2], tmp[3]] += a2 * b0
umlal	r11, r12, r5, r6
@ tmp[2] &= mask32
and.w	r11, r3, r11
@ tmp[3] &= mask32
and.w	r12, r3, r12
@ out[0] ^= (tmp[2]<<1)
eor.w	r2, r2, r11, lsl #1
@ 9, [tmp[3], tmp[4]] += a2 * b1
umlal	r12, r14, r5, r7
@ tmp[3] &= mask32
and.w	r12, r3, r12
@ 10, [tmp[4], tmp[5]] += a2 * b2
umlal	r14, r10, r5, r8
@ tmp[4] &= mask32
and.w	r14, r3, r14
@ a3 = (a >> 3) & mask32
and.w	r4, r3, r1, lsr #3
@ 11, [tmp[3], tmp[4]] += a3 * b0
umlal	r12, r14, r4, r6
@ tmp[3] &= mask32
and.w	r12, r3, r12
@ out[0] ^= (tmp[3]<<3)
eor.w	r2, r2, r12, lsl #3
@ 12, [tmp[6], tmp[7]] = a3 * b3
umull	r11, r12, r4, r9
@ tmp[7] &= mask32
and.w	r12, r3, r12
@ 13, [tmp[5], tmp[6]] += a2 * b3
umlal	r10, r11, r5, r9
@ tmp[5] &= mask32
and.w	r10, r3, r10
@ 14, [tmp[4], tmp[5]] += a3 * b1
umlal	r14, r10, r4, r7
@ tmp[4] &= mask32
and.w	r14, r3, r14
@ 15, [tmp[5], tmp[6]] += a3 * b2
umlal	r10, r11, r4, r8
@ tmp[5] &= mask32
and.w	r10, r3, r10
@ tmp[6] &= mask32
and.w	r11, r3, r11
@ out[1] = tmp[4] ^ (tmp[5]<<2)
eor.w	r3, r14, r10, lsl #2
@ out[1] ^= (tmp[6]<<1)
eor.w	r3, r3, r11, lsl #1
@ out[1] ^= (tmp[7]<<3)
eor.w	r3, r3, r12, lsl #3
@ store out[0]
str	r2, [r0, #0]
@ store out[1]
str	r3, [r0, #4]
ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}