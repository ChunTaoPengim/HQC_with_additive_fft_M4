.syntax unified
.cpu cortex-m4
.thumb
.fpu fpv4-sp-d16
.global radix_16_matrix_vector
.type radix_16_matrix_vector, %function
.align 2
radix_16_matrix_vector:

stmdb   sp!, { r4-r11, lr }
VSTMDB  sp!, { s0-s16 }

@ mconst = 0x11111111
mov.w	r11, #0x11111111
@ pconst = 0x00011101
movw    r12, #0x1101
movt    r12, #0x0001

@ r0: matrix, r1: vector

mov     r14, r0
vldm    r1, {s0-s15}

@ for i in range(16)
@ i = 0
ldr.w   r9, [r14, #0]
vmov    r10, s0
and.w   r8, r11, r9, lsr#0
umull   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umull   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umull   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umull   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 1
ldr.w   r9, [r14, #4]
vmov    r10, s1
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 2
ldr.w   r9, [r14, #8]
vmov    r10, s2
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 3
ldr.w   r9, [r14, #12]
vmov    r10, s3
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 4
ldr.w   r9, [r14, #16]
vmov    r10, s4
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 5
ldr.w   r9, [r14, #20]
vmov    r10, s5
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 6
ldr.w   r9, [r14, #24]
vmov    r10, s6
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 7
ldr.w   r9, [r14, #28]
vmov    r10, s7
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 8
ldr.w   r9, [r14, #32]
vmov    r10, s8
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 9
ldr.w   r9, [r14, #36]
vmov    r10, s9
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 10
ldr.w   r9, [r14, #40]
vmov    r10, s10
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 11
ldr.w   r9, [r14, #44]
vmov    r10, s11
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 12
ldr.w   r9, [r14, #48]
vmov    r10, s12
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 13
ldr.w   r9, [r14, #52]
vmov    r10, s13
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11
and.w   r7, r7, r11

@ i = 14
ldr.w   r9, [r14, #56]
vmov    r10, s14
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r4, r4, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r6, r6, r11

@ i = 15
ldr.w   r9, [r14, #60]
vmov    r10, s15
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r3, r3, r11
and.w   r8, r11, r9, lsr#1
umlal   r4, r5, r8, r10
and.w   r5, r5, r11
and.w   r8, r11, r9, lsr#3
umlal   r6, r7, r8, r10
and.w   r7, r7, r11

@ reduction
umull   r8, r9, r1, r12
mla     r8, r9, r12, r8
eor.w   r0, r0, r8
and.w   r0, r0, r11
umull   r8, r9, r3, r12
mla     r8, r9, r12, r8
eor.w   r2, r2, r8
and.w   r2, r2, r11
umull   r8, r9, r5, r12
mla     r8, r9, r12, r8
eor.w   r4, r4, r8
and.w   r4, r4, r11
umull   r8, r9, r7, r12
mla     r8, r9, r12, r8
eor.w   r6, r6, r8
and.w   r6, r6, r11

orr.w   r0, r0, r2, lsl#2
orr.w   r0, r0, r4, lsl#1
orr.w   r0, r0, r6, lsl#3

VLDMIA  sp!, { s0-s16 }
ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}