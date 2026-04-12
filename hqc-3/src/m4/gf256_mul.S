.syntax unified
.cpu cortex-m4
.thumb
.fpu fpv4-sp-d16
.global gf256_mul
.type gf256_mul, %function
.align 2
gf256_mul:
    push.w {r4-r11, r14}
    vpush.w {s16-s23}
    mov.w r12, #0x01010101
    mov.w r14, #0x1d

    @ precomputation of b' = b, bx, bx^2, ..., bx^7
    @ store them in s16 to s23

    @ bx^0
    vmov.w s16, r1
    @ bx^1
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s17, r1
    @ bx^2
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s18, r1
    @ bx^3
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s19, r1
    @ bx^4
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s20, r1
    @ bx^5
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s21, r1
    @ bx^6
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s22, r1
    @ bx^7
    and.w r3, r12, r1, lsr #7
    eor.w r1, r1, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r1, r3, r1, lsl#1
    vmov.w s23, r1

    @ gf256_mul

    @ load a0~a3 to r2~r5
    ldm r0, {r2-r5}

    @ assign acc0~acc3 to r6~r9
    @ tmp0~tmp1 to r10~r11

    @ k = 0
    vmov.w r11, s16
    and.w r10, r12, r2
    mul.w r6, r11, r10
    and.w r10, r12, r3
    mul.w r7, r11, r10
    and.w r10, r12, r4
    mul.w r8, r11, r10
    and.w r10, r12, r5
    mul.w r9, r11, r10
    
    @ k = 1
    vmov.w r11, s17
    and.w r10, r12, r2, lsr #1
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #1
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #1
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #1
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    @ k = 2
    vmov.w r11, s18
    and.w r10, r12, r2, lsr #2
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #2
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #2
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #2
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    @ k = 3
    vmov.w r11, s19
    and.w r10, r12, r2, lsr #3
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #3
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #3
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #3
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    @ k = 4
    vmov.w r11, s20
    and.w r10, r12, r2, lsr #4
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #4
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #4
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #4
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    @ k = 5
    vmov.w r11, s21
    and.w r10, r12, r2, lsr #5
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #5
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #5
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #5
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    @ k = 6
    vmov.w r11, s22
    and.w r10, r12, r2, lsr #6
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #6
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #6
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #6
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    @ k = 7
    vmov.w r11, s23
    and.w r10, r12, r2, lsr #7
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #7
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #7
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #7
    mul.w r10, r11, r10
    eor.w r9, r9, r10

    stm r0, {r6-r9}

    @ load a4~a7 to r2~r5
    @ load a8 to r14
    add r10, r0, #16
    ldm r10, {r2-r5, r14}

    @ assign acc0~acc3 to r6~r9
    @ tmp0~tmp1 to r10~r11

    @ k = 0
    vmov.w r11, s16
    and.w r10, r12, r2
    mul.w r6, r11, r10
    and.w r10, r12, r3
    mul.w r7, r11, r10
    and.w r10, r12, r4
    mul.w r8, r11, r10
    and.w r10, r12, r5
    mul.w r9, r11, r10
    and.w r10, r12, r14
    mul.w r1, r11, r10
    @ k = 1
    vmov.w r11, s17
    and.w r10, r12, r2, lsr #1
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #1
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #1
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #1
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #1
    mul.w r10, r11, r10
    eor.w r1, r1, r10
    @ k = 2
    vmov.w r11, s18
    and.w r10, r12, r2, lsr #2
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #2
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #2
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #2
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #2
    mul.w r10, r11, r10
    eor.w r1, r1, r10
    @ k = 3
    vmov.w r11, s19
    and.w r10, r12, r2, lsr #3
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #3
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #3
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #3
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #3
    mul.w r10, r11, r10
    eor.w r1, r1, r10
    @ k = 4
    vmov.w r11, s20
    and.w r10, r12, r2, lsr #4
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #4
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #4
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #4
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #4
    mul.w r10, r11, r10
    eor.w r1, r1, r10
    @ k = 5
    vmov.w r11, s21
    and.w r10, r12, r2, lsr #5
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #5
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #5
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #5
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #5
    mul.w r10, r11, r10
    eor.w r1, r1, r10
    @ k = 6
    vmov.w r11, s22
    and.w r10, r12, r2, lsr #6
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #6
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #6
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #6
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #6
    mul.w r10, r11, r10
    eor.w r1, r1, r10
    @ k = 7
    vmov.w r11, s23
    and.w r10, r12, r2, lsr #7
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3, lsr #7
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4, lsr #7
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5, lsr #7
    mul.w r10, r11, r10
    eor.w r9, r9, r10
    and.w r10, r12, r14, lsr #7
    mul.w r10, r11, r10
    eor.w r14, r1, r10

    add r10, r0, #16
    stm r10, {r6-r9, r14}

    vpop.w {s16-s23}
    pop {r4-r11, pc}
