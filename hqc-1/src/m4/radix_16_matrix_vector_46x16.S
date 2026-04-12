.syntax unified
.cpu cortex-m4
.thumb
.fpu fpv4-sp-d16
.global radix_16_matrix_vector_46x16
.type radix_16_matrix_vector_46x16, %function
.align 2
radix_16_matrix_vector_46x16:

stmdb   sp!, { r4-r11, lr }
VSTMDB  sp!, { s0-s16 }

@ mconst = 0x11111111
mov.w	r11, #0x11111111
@ pconst = 0x00011101
movw    r12, #0x1101
movt    r12, #0x0001

@ r0: input, r1: matrix, r2: vector

mov     r14, r1
vldm    r2, {s0-s15}

vmov    s16, r0

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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #0]

    
@ i = 0
ldr.w   r9, [r14, #64]
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
ldr.w   r9, [r14, #68]
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
ldr.w   r9, [r14, #72]
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
ldr.w   r9, [r14, #76]
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
ldr.w   r9, [r14, #80]
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
ldr.w   r9, [r14, #84]
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
ldr.w   r9, [r14, #88]
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
ldr.w   r9, [r14, #92]
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
ldr.w   r9, [r14, #96]
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
ldr.w   r9, [r14, #100]
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
ldr.w   r9, [r14, #104]
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
ldr.w   r9, [r14, #108]
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
ldr.w   r9, [r14, #112]
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
ldr.w   r9, [r14, #116]
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
ldr.w   r9, [r14, #120]
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
ldr.w   r9, [r14, #124]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #4]

    
@ i = 0
ldr.w   r9, [r14, #128]
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
ldr.w   r9, [r14, #132]
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
ldr.w   r9, [r14, #136]
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
ldr.w   r9, [r14, #140]
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
ldr.w   r9, [r14, #144]
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
ldr.w   r9, [r14, #148]
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
ldr.w   r9, [r14, #152]
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
ldr.w   r9, [r14, #156]
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
ldr.w   r9, [r14, #160]
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
ldr.w   r9, [r14, #164]
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
ldr.w   r9, [r14, #168]
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
ldr.w   r9, [r14, #172]
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
ldr.w   r9, [r14, #176]
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
ldr.w   r9, [r14, #180]
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
ldr.w   r9, [r14, #184]
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
ldr.w   r9, [r14, #188]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #8]

    
@ i = 0
ldr.w   r9, [r14, #192]
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
ldr.w   r9, [r14, #196]
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
ldr.w   r9, [r14, #200]
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
ldr.w   r9, [r14, #204]
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
ldr.w   r9, [r14, #208]
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
ldr.w   r9, [r14, #212]
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
ldr.w   r9, [r14, #216]
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
ldr.w   r9, [r14, #220]
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
ldr.w   r9, [r14, #224]
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
ldr.w   r9, [r14, #228]
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
ldr.w   r9, [r14, #232]
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
ldr.w   r9, [r14, #236]
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
ldr.w   r9, [r14, #240]
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
ldr.w   r9, [r14, #244]
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
ldr.w   r9, [r14, #248]
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
ldr.w   r9, [r14, #252]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #12]

    
@ i = 0
ldr.w   r9, [r14, #256]
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
ldr.w   r9, [r14, #260]
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
ldr.w   r9, [r14, #264]
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
ldr.w   r9, [r14, #268]
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
ldr.w   r9, [r14, #272]
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
ldr.w   r9, [r14, #276]
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
ldr.w   r9, [r14, #280]
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
ldr.w   r9, [r14, #284]
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
ldr.w   r9, [r14, #288]
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
ldr.w   r9, [r14, #292]
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
ldr.w   r9, [r14, #296]
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
ldr.w   r9, [r14, #300]
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
ldr.w   r9, [r14, #304]
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
ldr.w   r9, [r14, #308]
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
ldr.w   r9, [r14, #312]
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
ldr.w   r9, [r14, #316]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #16]

    
@ i = 0
ldr.w   r9, [r14, #320]
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
ldr.w   r9, [r14, #324]
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
ldr.w   r9, [r14, #328]
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
ldr.w   r9, [r14, #332]
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
ldr.w   r9, [r14, #336]
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
ldr.w   r9, [r14, #340]
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
ldr.w   r9, [r14, #344]
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
ldr.w   r9, [r14, #348]
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
ldr.w   r9, [r14, #352]
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
ldr.w   r9, [r14, #356]
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
ldr.w   r9, [r14, #360]
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
ldr.w   r9, [r14, #364]
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
ldr.w   r9, [r14, #368]
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
ldr.w   r9, [r14, #372]
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
ldr.w   r9, [r14, #376]
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
ldr.w   r9, [r14, #380]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #20]

    
@ i = 0
ldr.w   r9, [r14, #384]
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
ldr.w   r9, [r14, #388]
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
ldr.w   r9, [r14, #392]
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
ldr.w   r9, [r14, #396]
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
ldr.w   r9, [r14, #400]
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
ldr.w   r9, [r14, #404]
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
ldr.w   r9, [r14, #408]
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
ldr.w   r9, [r14, #412]
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
ldr.w   r9, [r14, #416]
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
ldr.w   r9, [r14, #420]
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
ldr.w   r9, [r14, #424]
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
ldr.w   r9, [r14, #428]
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
ldr.w   r9, [r14, #432]
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
ldr.w   r9, [r14, #436]
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
ldr.w   r9, [r14, #440]
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
ldr.w   r9, [r14, #444]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #24]

    
@ i = 0
ldr.w   r9, [r14, #448]
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
ldr.w   r9, [r14, #452]
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
ldr.w   r9, [r14, #456]
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
ldr.w   r9, [r14, #460]
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
ldr.w   r9, [r14, #464]
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
ldr.w   r9, [r14, #468]
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
ldr.w   r9, [r14, #472]
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
ldr.w   r9, [r14, #476]
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
ldr.w   r9, [r14, #480]
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
ldr.w   r9, [r14, #484]
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
ldr.w   r9, [r14, #488]
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
ldr.w   r9, [r14, #492]
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
ldr.w   r9, [r14, #496]
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
ldr.w   r9, [r14, #500]
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
ldr.w   r9, [r14, #504]
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
ldr.w   r9, [r14, #508]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #28]

    
@ i = 0
ldr.w   r9, [r14, #512]
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
ldr.w   r9, [r14, #516]
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
ldr.w   r9, [r14, #520]
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
ldr.w   r9, [r14, #524]
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
ldr.w   r9, [r14, #528]
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
ldr.w   r9, [r14, #532]
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
ldr.w   r9, [r14, #536]
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
ldr.w   r9, [r14, #540]
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
ldr.w   r9, [r14, #544]
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
ldr.w   r9, [r14, #548]
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
ldr.w   r9, [r14, #552]
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
ldr.w   r9, [r14, #556]
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
ldr.w   r9, [r14, #560]
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
ldr.w   r9, [r14, #564]
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
ldr.w   r9, [r14, #568]
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
ldr.w   r9, [r14, #572]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #32]

    
@ i = 0
ldr.w   r9, [r14, #576]
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
ldr.w   r9, [r14, #580]
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
ldr.w   r9, [r14, #584]
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
ldr.w   r9, [r14, #588]
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
ldr.w   r9, [r14, #592]
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
ldr.w   r9, [r14, #596]
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
ldr.w   r9, [r14, #600]
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
ldr.w   r9, [r14, #604]
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
ldr.w   r9, [r14, #608]
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
ldr.w   r9, [r14, #612]
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
ldr.w   r9, [r14, #616]
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
ldr.w   r9, [r14, #620]
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
ldr.w   r9, [r14, #624]
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
ldr.w   r9, [r14, #628]
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
ldr.w   r9, [r14, #632]
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
ldr.w   r9, [r14, #636]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #36]

    
@ i = 0
ldr.w   r9, [r14, #640]
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
ldr.w   r9, [r14, #644]
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
ldr.w   r9, [r14, #648]
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
ldr.w   r9, [r14, #652]
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
ldr.w   r9, [r14, #656]
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
ldr.w   r9, [r14, #660]
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
ldr.w   r9, [r14, #664]
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
ldr.w   r9, [r14, #668]
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
ldr.w   r9, [r14, #672]
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
ldr.w   r9, [r14, #676]
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
ldr.w   r9, [r14, #680]
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
ldr.w   r9, [r14, #684]
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
ldr.w   r9, [r14, #688]
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
ldr.w   r9, [r14, #692]
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
ldr.w   r9, [r14, #696]
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
ldr.w   r9, [r14, #700]
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

orr.w   r1, r0, r2, lsl#2
orr.w   r1, r1, r4, lsl#1
orr.w   r1, r1, r6, lsl#3

vmov    r0, s16
str     r1, [r0, #40]

    
ldr.w   r9, [r14, #704]
vmov    r10, s0
and.w   r8, r11, r9, lsr#0
umull   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umull   r2, r3, r8, r10
and.w   r2, r2, r11

@ i = 1
ldr.w   r9, [r14, #708]
vmov    r10, s1
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 2
ldr.w   r9, [r14, #712]
vmov    r10, s2
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 3
ldr.w   r9, [r14, #716]
vmov    r10, s3
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 4
ldr.w   r9, [r14, #720]
vmov    r10, s4
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 5
ldr.w   r9, [r14, #724]
vmov    r10, s5
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 6
ldr.w   r9, [r14, #728]
vmov    r10, s6
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 7
ldr.w   r9, [r14, #732]
vmov    r10, s7
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 8
ldr.w   r9, [r14, #736]
vmov    r10, s8
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 9
ldr.w   r9, [r14, #740]
vmov    r10, s9
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 10
ldr.w   r9, [r14, #744]
vmov    r10, s10
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 11
ldr.w   r9, [r14, #748]
vmov    r10, s11
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 12
ldr.w   r9, [r14, #752]
vmov    r10, s12
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 13
ldr.w   r9, [r14, #756]
vmov    r10, s13
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11
and.w   r3, r3, r11


@ i = 14
ldr.w   r9, [r14, #760]
vmov    r10, s14
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r0, r0, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r2, r2, r11


@ i = 15
ldr.w   r9, [r14, #764]
vmov    r10, s15
and.w   r8, r11, r9, lsr#0
umlal   r0, r1, r8, r10
and.w   r1, r1, r11
and.w   r8, r11, r9, lsr#2
umlal   r2, r3, r8, r10
and.w   r3, r3, r11

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

orr.w   r1, r0, r2, lsl#2
          
vmov    r0, s16
str     r1, [r0, #44]

    
VLDMIA  sp!, { s0-s16 }
ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}

