    .syntax unified
    .cpu cortex-m4
    .thumb
    .fpu fpv4-sp-d16
    .global radix_16_matrix_vector_32x56
    .type radix_16_matrix_vector_32x56, %function
    .align 2
    radix_16_matrix_vector_32x56:

    stmdb   sp!, { r4-r11, lr }
    VSTMDB  sp!, { s0-s29 }

    @ mconst = 0x11111111
    mov.w	r11, #0x11111111
    @ pconst = 0x00011101
    movw    r12, #0x1101
    movt    r12, #0x0001

    @ r0: input, r1: matrix, r2: vector
    mov     r14, r1
    vmov    s28, r0
    
    vldm    r2, {s0-s27}
    vmov    s29, r2
            
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

    @ i = 16
    ldr.w   r9, [r14, #64]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #68]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #72]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #76]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #80]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #84]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #88]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #92]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #96]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #100]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #104]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #108]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #0]
                
    @ i = 0
    ldr.w   r9, [r14, #224]
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
    ldr.w   r9, [r14, #228]
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
    ldr.w   r9, [r14, #232]
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
    ldr.w   r9, [r14, #236]
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
    ldr.w   r9, [r14, #240]
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
    ldr.w   r9, [r14, #244]
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
    ldr.w   r9, [r14, #248]
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
    ldr.w   r9, [r14, #252]
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
    ldr.w   r9, [r14, #256]
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
    ldr.w   r9, [r14, #260]
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
    ldr.w   r9, [r14, #264]
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
    ldr.w   r9, [r14, #268]
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
    ldr.w   r9, [r14, #272]
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
    ldr.w   r9, [r14, #276]
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
    ldr.w   r9, [r14, #280]
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
    ldr.w   r9, [r14, #284]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #288]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #292]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #296]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #300]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #304]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #308]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #312]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #316]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #320]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #324]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #328]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #332]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #4]
                
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

    @ i = 16
    ldr.w   r9, [r14, #512]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #516]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #520]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #524]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #528]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #532]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #536]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #540]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #544]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #548]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #552]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #556]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #8]
                
    @ i = 0
    ldr.w   r9, [r14, #672]
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
    ldr.w   r9, [r14, #676]
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
    ldr.w   r9, [r14, #680]
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
    ldr.w   r9, [r14, #684]
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
    ldr.w   r9, [r14, #688]
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
    ldr.w   r9, [r14, #692]
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
    ldr.w   r9, [r14, #696]
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
    ldr.w   r9, [r14, #700]
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
    ldr.w   r9, [r14, #704]
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
    ldr.w   r9, [r14, #708]
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
    ldr.w   r9, [r14, #712]
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
    ldr.w   r9, [r14, #716]
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
    ldr.w   r9, [r14, #720]
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
    ldr.w   r9, [r14, #724]
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
    ldr.w   r9, [r14, #728]
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
    ldr.w   r9, [r14, #732]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #736]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #740]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #744]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #748]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #752]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #756]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #760]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #764]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #768]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #772]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #776]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #780]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #12]
                
    @ i = 0
    ldr.w   r9, [r14, #896]
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
    ldr.w   r9, [r14, #900]
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
    ldr.w   r9, [r14, #904]
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
    ldr.w   r9, [r14, #908]
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
    ldr.w   r9, [r14, #912]
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
    ldr.w   r9, [r14, #916]
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
    ldr.w   r9, [r14, #920]
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
    ldr.w   r9, [r14, #924]
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
    ldr.w   r9, [r14, #928]
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
    ldr.w   r9, [r14, #932]
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
    ldr.w   r9, [r14, #936]
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
    ldr.w   r9, [r14, #940]
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
    ldr.w   r9, [r14, #944]
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
    ldr.w   r9, [r14, #948]
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
    ldr.w   r9, [r14, #952]
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
    ldr.w   r9, [r14, #956]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #960]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #964]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #968]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #972]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #976]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #980]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #984]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #988]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #992]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #996]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1000]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1004]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #16]
                
    @ i = 0
    ldr.w   r9, [r14, #1120]
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
    ldr.w   r9, [r14, #1124]
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
    ldr.w   r9, [r14, #1128]
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
    ldr.w   r9, [r14, #1132]
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
    ldr.w   r9, [r14, #1136]
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
    ldr.w   r9, [r14, #1140]
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
    ldr.w   r9, [r14, #1144]
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
    ldr.w   r9, [r14, #1148]
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
    ldr.w   r9, [r14, #1152]
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
    ldr.w   r9, [r14, #1156]
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
    ldr.w   r9, [r14, #1160]
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
    ldr.w   r9, [r14, #1164]
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
    ldr.w   r9, [r14, #1168]
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
    ldr.w   r9, [r14, #1172]
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
    ldr.w   r9, [r14, #1176]
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
    ldr.w   r9, [r14, #1180]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1184]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1188]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1192]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1196]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1200]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1204]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1208]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1212]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1216]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1220]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1224]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1228]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #20]
                
    @ i = 0
    ldr.w   r9, [r14, #1344]
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
    ldr.w   r9, [r14, #1348]
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
    ldr.w   r9, [r14, #1352]
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
    ldr.w   r9, [r14, #1356]
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
    ldr.w   r9, [r14, #1360]
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
    ldr.w   r9, [r14, #1364]
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
    ldr.w   r9, [r14, #1368]
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
    ldr.w   r9, [r14, #1372]
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
    ldr.w   r9, [r14, #1376]
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
    ldr.w   r9, [r14, #1380]
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
    ldr.w   r9, [r14, #1384]
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
    ldr.w   r9, [r14, #1388]
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
    ldr.w   r9, [r14, #1392]
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
    ldr.w   r9, [r14, #1396]
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
    ldr.w   r9, [r14, #1400]
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
    ldr.w   r9, [r14, #1404]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1408]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1412]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1416]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1420]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1424]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1428]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1432]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1436]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1440]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1444]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1448]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1452]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #24]
                
    @ i = 0
    ldr.w   r9, [r14, #1568]
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
    ldr.w   r9, [r14, #1572]
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
    ldr.w   r9, [r14, #1576]
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
    ldr.w   r9, [r14, #1580]
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
    ldr.w   r9, [r14, #1584]
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
    ldr.w   r9, [r14, #1588]
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
    ldr.w   r9, [r14, #1592]
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
    ldr.w   r9, [r14, #1596]
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
    ldr.w   r9, [r14, #1600]
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
    ldr.w   r9, [r14, #1604]
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
    ldr.w   r9, [r14, #1608]
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
    ldr.w   r9, [r14, #1612]
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
    ldr.w   r9, [r14, #1616]
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
    ldr.w   r9, [r14, #1620]
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
    ldr.w   r9, [r14, #1624]
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
    ldr.w   r9, [r14, #1628]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1632]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1636]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1640]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1644]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1648]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1652]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1656]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1660]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1664]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1668]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1672]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1676]
    vmov    r10, s27
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
        
    vmov    r0, s28
    str     r1, [r0, #28]
                
    vmov    r2, s29
    add     r2, #112
    vldm    r2, {s0-s27}
    vmov    s29, r2
            
    @ i = 0
    ldr.w   r9, [r14, #112]
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
    ldr.w   r9, [r14, #116]
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
    ldr.w   r9, [r14, #120]
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
    ldr.w   r9, [r14, #124]
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
    ldr.w   r9, [r14, #128]
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
    ldr.w   r9, [r14, #132]
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
    ldr.w   r9, [r14, #136]
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
    ldr.w   r9, [r14, #140]
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
    ldr.w   r9, [r14, #144]
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
    ldr.w   r9, [r14, #148]
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
    ldr.w   r9, [r14, #152]
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
    ldr.w   r9, [r14, #156]
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
    ldr.w   r9, [r14, #160]
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
    ldr.w   r9, [r14, #164]
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
    ldr.w   r9, [r14, #168]
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
    ldr.w   r9, [r14, #172]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #176]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #180]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #184]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #188]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #192]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #196]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #200]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #204]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #208]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #212]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #216]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #220]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #0]
    eor.w   r1, r1, r2
    str     r1, [r0, #0]
                
    @ i = 0
    ldr.w   r9, [r14, #336]
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
    ldr.w   r9, [r14, #340]
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
    ldr.w   r9, [r14, #344]
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
    ldr.w   r9, [r14, #348]
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
    ldr.w   r9, [r14, #352]
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
    ldr.w   r9, [r14, #356]
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
    ldr.w   r9, [r14, #360]
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
    ldr.w   r9, [r14, #364]
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
    ldr.w   r9, [r14, #368]
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
    ldr.w   r9, [r14, #372]
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
    ldr.w   r9, [r14, #376]
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
    ldr.w   r9, [r14, #380]
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
    ldr.w   r9, [r14, #384]
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
    ldr.w   r9, [r14, #388]
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
    ldr.w   r9, [r14, #392]
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
    ldr.w   r9, [r14, #396]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #400]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #404]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #408]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #412]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #416]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #420]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #424]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #428]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #432]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #436]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #440]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #444]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #4]
    eor.w   r1, r1, r2
    str     r1, [r0, #4]
                
    @ i = 0
    ldr.w   r9, [r14, #560]
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
    ldr.w   r9, [r14, #564]
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
    ldr.w   r9, [r14, #568]
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
    ldr.w   r9, [r14, #572]
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
    ldr.w   r9, [r14, #576]
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
    ldr.w   r9, [r14, #580]
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
    ldr.w   r9, [r14, #584]
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
    ldr.w   r9, [r14, #588]
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
    ldr.w   r9, [r14, #592]
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
    ldr.w   r9, [r14, #596]
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
    ldr.w   r9, [r14, #600]
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
    ldr.w   r9, [r14, #604]
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
    ldr.w   r9, [r14, #608]
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
    ldr.w   r9, [r14, #612]
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
    ldr.w   r9, [r14, #616]
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
    ldr.w   r9, [r14, #620]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #624]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #628]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #632]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #636]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #640]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #644]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #648]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #652]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #656]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #660]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #664]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #668]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #8]
    eor.w   r1, r1, r2
    str     r1, [r0, #8]
                
    @ i = 0
    ldr.w   r9, [r14, #784]
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
    ldr.w   r9, [r14, #788]
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
    ldr.w   r9, [r14, #792]
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
    ldr.w   r9, [r14, #796]
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
    ldr.w   r9, [r14, #800]
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
    ldr.w   r9, [r14, #804]
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
    ldr.w   r9, [r14, #808]
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
    ldr.w   r9, [r14, #812]
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
    ldr.w   r9, [r14, #816]
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
    ldr.w   r9, [r14, #820]
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
    ldr.w   r9, [r14, #824]
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
    ldr.w   r9, [r14, #828]
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
    ldr.w   r9, [r14, #832]
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
    ldr.w   r9, [r14, #836]
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
    ldr.w   r9, [r14, #840]
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
    ldr.w   r9, [r14, #844]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #848]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #852]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #856]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #860]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #864]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #868]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #872]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #876]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #880]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #884]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #888]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #892]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #12]
    eor.w   r1, r1, r2
    str     r1, [r0, #12]
                
    @ i = 0
    ldr.w   r9, [r14, #1008]
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
    ldr.w   r9, [r14, #1012]
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
    ldr.w   r9, [r14, #1016]
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
    ldr.w   r9, [r14, #1020]
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
    ldr.w   r9, [r14, #1024]
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
    ldr.w   r9, [r14, #1028]
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
    ldr.w   r9, [r14, #1032]
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
    ldr.w   r9, [r14, #1036]
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
    ldr.w   r9, [r14, #1040]
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
    ldr.w   r9, [r14, #1044]
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
    ldr.w   r9, [r14, #1048]
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
    ldr.w   r9, [r14, #1052]
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
    ldr.w   r9, [r14, #1056]
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
    ldr.w   r9, [r14, #1060]
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
    ldr.w   r9, [r14, #1064]
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
    ldr.w   r9, [r14, #1068]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1072]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1076]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1080]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1084]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1088]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1092]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1096]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1100]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1104]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1108]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1112]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1116]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #16]
    eor.w   r1, r1, r2
    str     r1, [r0, #16]
                
    @ i = 0
    ldr.w   r9, [r14, #1232]
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
    ldr.w   r9, [r14, #1236]
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
    ldr.w   r9, [r14, #1240]
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
    ldr.w   r9, [r14, #1244]
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
    ldr.w   r9, [r14, #1248]
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
    ldr.w   r9, [r14, #1252]
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
    ldr.w   r9, [r14, #1256]
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
    ldr.w   r9, [r14, #1260]
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
    ldr.w   r9, [r14, #1264]
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
    ldr.w   r9, [r14, #1268]
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
    ldr.w   r9, [r14, #1272]
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
    ldr.w   r9, [r14, #1276]
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
    ldr.w   r9, [r14, #1280]
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
    ldr.w   r9, [r14, #1284]
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
    ldr.w   r9, [r14, #1288]
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
    ldr.w   r9, [r14, #1292]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1296]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1300]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1304]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1308]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1312]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1316]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1320]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1324]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1328]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1332]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1336]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1340]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #20]
    eor.w   r1, r1, r2
    str     r1, [r0, #20]
                
    @ i = 0
    ldr.w   r9, [r14, #1456]
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
    ldr.w   r9, [r14, #1460]
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
    ldr.w   r9, [r14, #1464]
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
    ldr.w   r9, [r14, #1468]
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
    ldr.w   r9, [r14, #1472]
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
    ldr.w   r9, [r14, #1476]
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
    ldr.w   r9, [r14, #1480]
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
    ldr.w   r9, [r14, #1484]
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
    ldr.w   r9, [r14, #1488]
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
    ldr.w   r9, [r14, #1492]
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
    ldr.w   r9, [r14, #1496]
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
    ldr.w   r9, [r14, #1500]
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
    ldr.w   r9, [r14, #1504]
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
    ldr.w   r9, [r14, #1508]
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
    ldr.w   r9, [r14, #1512]
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
    ldr.w   r9, [r14, #1516]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1520]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1524]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1528]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1532]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1536]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1540]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1544]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1548]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1552]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1556]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1560]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1564]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #24]
    eor.w   r1, r1, r2
    str     r1, [r0, #24]
                
    @ i = 0
    ldr.w   r9, [r14, #1680]
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
    ldr.w   r9, [r14, #1684]
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
    ldr.w   r9, [r14, #1688]
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
    ldr.w   r9, [r14, #1692]
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
    ldr.w   r9, [r14, #1696]
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
    ldr.w   r9, [r14, #1700]
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
    ldr.w   r9, [r14, #1704]
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
    ldr.w   r9, [r14, #1708]
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
    ldr.w   r9, [r14, #1712]
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
    ldr.w   r9, [r14, #1716]
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
    ldr.w   r9, [r14, #1720]
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
    ldr.w   r9, [r14, #1724]
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
    ldr.w   r9, [r14, #1728]
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
    ldr.w   r9, [r14, #1732]
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
    ldr.w   r9, [r14, #1736]
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
    ldr.w   r9, [r14, #1740]
    vmov    r10, s15
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

    @ i = 16
    ldr.w   r9, [r14, #1744]
    vmov    r10, s16
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

    @ i = 17
    ldr.w   r9, [r14, #1748]
    vmov    r10, s17
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

    @ i = 18
    ldr.w   r9, [r14, #1752]
    vmov    r10, s18
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

    @ i = 19
    ldr.w   r9, [r14, #1756]
    vmov    r10, s19
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

    @ i = 20
    ldr.w   r9, [r14, #1760]
    vmov    r10, s20
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

    @ i = 21
    ldr.w   r9, [r14, #1764]
    vmov    r10, s21
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

    @ i = 22
    ldr.w   r9, [r14, #1768]
    vmov    r10, s22
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

    @ i = 23
    ldr.w   r9, [r14, #1772]
    vmov    r10, s23
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

    @ i = 24
    ldr.w   r9, [r14, #1776]
    vmov    r10, s24
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
    
    @ i = 25
    ldr.w   r9, [r14, #1780]
    vmov    r10, s25
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

    @ i = 26
    ldr.w   r9, [r14, #1784]
    vmov    r10, s26
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
    
    @ i = 27
    ldr.w   r9, [r14, #1788]
    vmov    r10, s27
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
        
    vmov    r0, s28
    ldr     r2, [r0, #28]
    eor.w   r1, r1, r2
    str     r1, [r0, #28]
                
    VLDMIA  sp!, { s0-s29 }
    ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}
    
