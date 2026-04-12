// implementaed by daniel
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global mul2_asm
.type mul2_asm, %function
.align 2
mul2_asm:

stmdb   sp!, { r4-r11, lr }
VSTMDB  sp!, {s16-s21}

@ mask32 = 0x11111111
mov.w	r3, #286331153

@ store a and b in FP registers
@ a[1] -> s16
ldr r12, [r1, #4]
vmov s16, r12
@ b[1] -> s17
ldr r14, [r2, #4]
vmov s17, r14

@ a[0]
ldr r1, [r1, #0]
@ b[0]
ldr r2, [r2, #0]

@ store hs0 and hs1 in FP registers
@ hs0 -> s18
eor.w r12, r1, r12
vmov s18, r12
@ hs1 -> s19
eor.w r14, r2, r14
vmov s19, r14

@ a0 = (a[0])&mask32;
and.w	r4, r3, r1
@ b0 = (b[0])&mask32;
and.w	r6, r3, r2
@ b1 = (b[0] >> 2) & mask32
and.w	r7, r3, r2, lsr #2
@ b2 = (b[0] >> 1) & mask32
and.w	r8, r3, r2, lsr #1
@ b3 = (b[0] >> 3) & mask32
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
@ a1 = (a[0] >> 2) & mask32
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
@ a2 = (a[0] >> 1) & mask32
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
@ a3 = (a[0] >> 3) & mask32
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
eor.w	r4, r14, r10, lsl #2
@ out[1] ^= (tmp[6]<<1)
eor.w	r4, r4, r11, lsl #1
@ out[1] ^= (tmp[7]<<3)
eor.w	r4, r4, r12, lsl #3
@ store out[0] in s20
@ str	r2, [r0, #0]
vmov s20, r2
@ store out[1] in s21
@ str	r4, [r0, #4]
vmov s21, r4

@ move back a and b from FP registers
vmov r1, s16
vmov r2, s17
@ a0 = (a[1])&mask32;
and.w	r4, r3, r1
@ b0 = (b[1])&mask32;
and.w	r6, r3, r2
@ b1 = (b[1] >> 2) & mask32
and.w	r7, r3, r2, lsr #2
@ b2 = (b[1] >> 1) & mask32
and.w	r8, r3, r2, lsr #1
@ b3 = (b[1] >> 3) & mask32
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
@ a1 = (a[1] >> 2) & mask32
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
@ a2 = (a[1] >> 1) & mask32
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
@ a3 = (a[1] >> 3) & mask32
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
eor.w	r4, r14, r10, lsl #2
@ out[1] ^= (tmp[6]<<1)
eor.w	r4, r4, r11, lsl #1
@ out[1] ^= (tmp[7]<<3)
eor.w	r4, r4, r12, lsl #3

@ move back hs0 and hs1 from FP registers; store out[0] in s18 and store out[1] in s19
vmov r1, s18
vmov s18, r2
vmov r2, s19
vmov s19, r4
@ a0 = (hs0)&mask32;
and.w	r4, r3, r1
@ b0 = (hs1)&mask32;
and.w	r6, r3, r2
@ b1 = (hs1 >> 2) & mask32
and.w	r7, r3, r2, lsr #2
@ b2 = (hs1 >> 1) & mask32
and.w	r8, r3, r2, lsr #1
@ b3 = (hs1 >> 3) & mask32
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
@ a1 = (hs0 >> 2) & mask32
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
@ a2 = (hs0 >> 1) & mask32
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
@ a3 = (hs0 >> 3) & mask32
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
eor.w	r4, r14, r10, lsl #2
@ out[1] ^= (tmp[6]<<1)
eor.w	r4, r4, r11, lsl #1
@ out[1] ^= (tmp[7]<<3)
eor.w	r4, r4, r12, lsl #3
@ here r2 is hl2[0], r4 is hl2[1]

@ move back c[0] to c[3]
vmov r5, s20
vmov r6, s21
vmov r7, s18
vmov r8, s19
@ hl2[0] = hl2[0] ^ c[0] ^ c[2]
eor.w r2, r2, r5
eor.w r2, r2, r7
@ hl2[1] = hl2[1] ^ c[1] ^ c[3]
eor.w r4, r4, r6
eor.w r4, r4, r8
@ c[1] ^= hl2[0]
eor.w r6, r2, r6
@ c[2] ^= hl2[1]
eor.w r7, r4, r7

stm r0, {r5-r8}

VLDMIA  sp!, {s16-s21}
ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}