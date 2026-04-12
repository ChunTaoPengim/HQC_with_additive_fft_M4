.syntax unified
.cpu cortex-m4
.thumb
.fpu fpv4-sp-d16
.global gf256_madd
.type gf256_madd, %function
.align 2
gf256_madd:
    push.w {r4-r11, r14}
    vpush.w {s16-s23}
    mov.w r12, #0x01010101
    mov.w r14, #0x1d

    @ precomputation of b' = b, bx, bx^2, ..., bx^7
    @ store them in s16 to s23

    @ bx^0
    vmov.w s16, r2
    @ bx^1
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s17, r2
    @ bx^2
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s18, r2
    @ bx^3
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s19, r2
    @ bx^4
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s20, r2
    @ bx^5
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s21, r2
    @ bx^6
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s22, r2
    @ bx^7
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s23, r2

    @ store r1 in s24
    vmov.w s24, r1

    @ gf256_madd

    @ load a0~a3 to r2~r5
    @ load a4 to r14
    ldm r1, {r2-r5, r14}

    @ assign acc0~acc3 to r6~r9
    @ assign acc4 to r1
    ldm r0, {r6-r10}
    
    @ tmp0,tmp1 to r1,r11

    @ k = 0
    vmov.w r11, s16
    and.w r1, r12, r2
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 1
    vmov.w r11, s17
    and.w r1, r12, r2, lsr #1
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #1
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #1
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #1
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #1
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 2
    vmov.w r11, s18
    and.w r1, r12, r2, lsr #2
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #2
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #2
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #2
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #2
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 3
    vmov.w r11, s19
    and.w r1, r12, r2, lsr #3
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #3
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #3
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #3
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #3
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 4
    vmov.w r11, s20
    and.w r1, r12, r2, lsr #4
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #4
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #4
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #4
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #4
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 5
    vmov.w r11, s21
    and.w r1, r12, r2, lsr #5
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #5
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #5
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #5
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #5
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 6
    vmov.w r11, s22
    and.w r1, r12, r2, lsr #6
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #6
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #6
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #6
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #6
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 7
    vmov.w r11, s23
    and.w r1, r12, r2, lsr #7
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #7
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #7
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #7
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #7
    mul.w r1, r11, r1
    eor.w r10, r10, r1

    stm r0, {r6-r10}

    @ load a5~a8 to r2~r5
    @ load a9 to r14
    vmov.w r1, s24
    add r10, r1, #20
    ldm r10, {r2-r5, r14}

    @ assign acc5~acc8 to r6~r9
    @ assign acc9 to r1
    add r11, r0, #20
    ldm r11, {r6-r10}
    
    @ tmp0,tmp1 to r1,r11

    @ k = 0
    vmov.w r11, s16
    and.w r1, r12, r2
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 1
    vmov.w r11, s17
    and.w r1, r12, r2, lsr #1
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #1
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #1
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #1
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #1
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 2
    vmov.w r11, s18
    and.w r1, r12, r2, lsr #2
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #2
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #2
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #2
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #2
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 3
    vmov.w r11, s19
    and.w r1, r12, r2, lsr #3
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #3
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #3
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #3
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #3
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 4
    vmov.w r11, s20
    and.w r1, r12, r2, lsr #4
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #4
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #4
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #4
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #4
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 5
    vmov.w r11, s21
    and.w r1, r12, r2, lsr #5
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #5
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #5
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #5
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #5
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 6
    vmov.w r11, s22
    and.w r1, r12, r2, lsr #6
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #6
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #6
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #6
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #6
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 7
    vmov.w r11, s23
    and.w r1, r12, r2, lsr #7
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #7
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #7
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #7
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #7
    mul.w r1, r11, r1
    eor.w r10, r10, r1

    add r11, r0, #20
    stm r11, {r6-r10}

    @ load a10~a13 to r2~r5
    @ load a14 to r14
    vmov.w r1, s24
    add r10, r1, #40
    ldm r10, {r2-r5, r14}

    @ assign acc10~acc13 to r6~r9
    @ assign acc14 to r1
    add r11, r0, #40
    ldm r11, {r6-r10}
    
    @ tmp0,tmp1 to r1,r11

    @ k = 0
    vmov.w r11, s16
    and.w r1, r12, r2
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 1
    vmov.w r11, s17
    and.w r1, r12, r2, lsr #1
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #1
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #1
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #1
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #1
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 2
    vmov.w r11, s18
    and.w r1, r12, r2, lsr #2
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #2
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #2
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #2
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #2
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 3
    vmov.w r11, s19
    and.w r1, r12, r2, lsr #3
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #3
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #3
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #3
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #3
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 4
    vmov.w r11, s20
    and.w r1, r12, r2, lsr #4
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #4
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #4
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #4
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #4
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 5
    vmov.w r11, s21
    and.w r1, r12, r2, lsr #5
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #5
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #5
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #5
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #5
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 6
    vmov.w r11, s22
    and.w r1, r12, r2, lsr #6
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #6
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #6
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #6
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #6
    mul.w r1, r11, r1
    eor.w r10, r10, r1
    @ k = 7
    vmov.w r11, s23
    and.w r1, r12, r2, lsr #7
    mul.w r1, r11, r1
    eor.w r6, r6, r1
    and.w r1, r12, r3, lsr #7
    mul.w r1, r11, r1
    eor.w r7, r7, r1
    and.w r1, r12, r4, lsr #7
    mul.w r1, r11, r1
    eor.w r8, r8, r1
    and.w r1, r12, r5, lsr #7
    mul.w r1, r11, r1
    eor.w r9, r9, r1
    and.w r1, r12, r14, lsr #7
    mul.w r1, r11, r1
    eor.w r10, r10, r1

    add r11, r0, #40
    stm r11, {r6-r10}

    vpop.w {s16-s23}
    pop {r4-r11, pc}

.syntax unified
.cpu cortex-m4
.thumb
.fpu fpv4-sp-d16
.global gf256_madd_32
.type gf256_madd, %function
.align 2
gf256_madd_32:
    push.w {r4-r11, r14}
    vpush.w {s16-s23}
    mov.w r12, #0x01010101
    mov.w r14, #0x1d

    @ precomputation of b' = b, bx, bx^2, ..., bx^7
    @ store them in s16 to s23

    @ bx^0
    vmov.w s16, r2
    @ bx^1
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s17, r2
    @ bx^2
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s18, r2
    @ bx^3
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s19, r2
    @ bx^4
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s20, r2
    @ bx^5
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s21, r2
    @ bx^6
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s22, r2
    @ bx^7
    and.w r3, r12, r2, lsr #7
    eor.w r2, r2, r3, lsl #7
    mul.w r3, r3, r14
    eor.w r2, r3, r2, lsl#1
    vmov.w s23, r2

    @ gf256_madd

    @ load a0~a3 to r2~r5
    ldm r1, {r2-r5}

    @ assign acc0~acc3 to r6~r9
    ldm r0, {r6-r9}
    
    @ tmp0~tmp1 to r10~r11

    @ k = 0
    vmov.w r11, s16
    and.w r10, r12, r2
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5
    mul.w r10, r11, r10
    eor.w r9, r9, r10
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
    add r10, r1, #16
    ldm r10, {r2-r5}

    @ assign acc4~acc7 to r6~r9
    add r10, r0, #16
    ldm r10, {r6-r9}
    
    @ tmp0~tmp1 to r10~r11

    @ k = 0
    vmov.w r11, s16
    and.w r10, r12, r2
    mul.w r10, r11, r10
    eor.w r6, r6, r10
    and.w r10, r12, r3
    mul.w r10, r11, r10
    eor.w r7, r7, r10
    and.w r10, r12, r4
    mul.w r10, r11, r10
    eor.w r8, r8, r10
    and.w r10, r12, r5
    mul.w r10, r11, r10
    eor.w r9, r9, r10
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

    add r10, r0, #16
    stm r10, {r6-r9}

    vpop.w {s16-s23}
    pop {r4-r11, pc}