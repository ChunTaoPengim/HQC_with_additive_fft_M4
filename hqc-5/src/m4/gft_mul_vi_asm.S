// Implemented by Myeonghoon Lee and Jihoon Jang
// public domain

.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v17
.type gft_mul_v17, %function
.align 2
gft_mul_v17:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldm r12, {s0-s31}
vmov r0, s2    @ load x[2] from s2
vmov r1, s5    @ load x[5] from s5
eor r0, r0, r1    @ x[2] = x[2] ^ x[5]
vmov r2, s25    @ load x[25] from s25
vmov r3, s29    @ load x[29] from s29
eor r2, r2, r3    @ x[25] = x[25] ^ x[29]
vmov r4, s27    @ load x[27] from s27
eor r4, r4, r0    @ x[27] = x[27] ^ x[2]
vmov r5, s31    @ load x[31] from s31
eor r4, r4, r5    @ x[27] = x[27] ^ x[31]
vmov r6, s7    @ load x[7] from s7
eor r6, r6, r2    @ x[7] = x[7] ^ x[25]
vmov r7, s3    @ load x[3] from s3
eor r7, r7, r4    @ x[3] = x[3] ^ x[27]
eor r5, r5, r6    @ x[31] = x[31] ^ x[7]
vmov r8, s22    @ load x[22] from s22
vmov r9, s18    @ load x[18] from s18
eor r8, r8, r9    @ x[22] = x[22] ^ x[18]
vmov r10, s20    @ load x[20] from s20
vmov r11, s28    @ load x[28] from s28
eor r10, r10, r11    @ x[20] = x[20] ^ x[28]
vmov s28, r11    @ spill x[28] from r11
vmov r11, s16    @ load x[16] from s16
eor r11, r11, r2    @ x[16] = x[16] ^ x[25]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s21    @ load x[21] from s21
eor r4, r4, r3    @ x[21] = x[21] ^ x[29]
vmov s29, r3    @ spill x[29] from r3
vmov r3, s11    @ load x[11] from s11
eor r0, r0, r3    @ x[2] = x[2] ^ x[11]
vmov s21, r4    @ spill x[21] from r4
vmov r4, s24    @ load x[24] from s24
eor r0, r0, r4    @ x[2] = x[2] ^ x[24]
vmov s2, r0    @ spill x[2] from r0
vmov r0, s9    @ load x[9] from s9
vmov s5, r1    @ spill x[5] from r1
vmov r1, s15    @ load x[15] from s15
eor r0, r0, r1    @ x[9] = x[9] ^ x[15]
eor r1, r1, r6    @ x[15] = x[15] ^ x[7]
eor r0, r0, r2    @ x[9] = x[9] ^ x[25]
vmov s15, r1    @ spill x[15] from r1
vmov r1, s26    @ load x[26] from s26
eor r1, r1, r7    @ x[26] = x[26] ^ x[3]
eor r4, r4, r6    @ x[24] = x[24] ^ x[7]
vmov s18, r9    @ spill x[18] from r9
vmov r9, s17    @ load x[17] from s17
eor r6, r6, r9    @ x[7] = x[7] ^ x[17]
eor r9, r9, r10    @ x[17] = x[17] ^ x[20]
eor r10, r10, r11    @ x[20] = x[20] ^ x[16]
vmov s17, r9    @ spill x[17] from r9
vmov r9, s0    @ load x[0] from s0
eor r9, r9, r3    @ x[0] = x[0] ^ x[11]
eor r9, r9, r5    @ x[0] = x[0] ^ x[31]
vmov s0, r9    @ spill x[0] from r9
vmov r9, s30    @ load x[30] from s30
eor r9, r9, r8    @ x[30] = x[30] ^ x[22]
eor r1, r1, r2    @ x[26] = x[26] ^ x[25]
vmov s26, r1    @ spill x[26] from r1
vmov r1, s1    @ load x[1] from s1
eor r9, r9, r1    @ x[30] = x[30] ^ x[1]
eor r0, r0, r7    @ x[9] = x[9] ^ x[3]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s6    @ load x[6] from s6
eor r9, r9, r7    @ x[6] = x[6] ^ x[3]
vmov s9, r0    @ spill x[9] from r0
vmov r0, s14    @ load x[14] from s14
eor r0, r0, r1    @ x[14] = x[14] ^ x[1]
eor r4, r4, r7    @ x[24] = x[24] ^ x[3]
eor r7, r7, r5    @ x[3] = x[3] ^ x[31]
vmov s3, r7    @ spill x[3] from r7
vmov r7, s23    @ load x[23] from s23
eor r5, r5, r7    @ x[31] = x[31] ^ x[23]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s18    @ load x[18] from s18
eor r5, r5, r4    @ x[31] = x[31] ^ x[18]
vmov s6, r9    @ spill x[6] from r9
vmov r9, s5    @ load x[5] from s5
eor r4, r4, r9    @ x[18] = x[18] ^ x[5]
vmov s18, r4    @ spill x[18] from r4
vmov r4, s21    @ load x[21] from s21
eor r9, r9, r4    @ x[5] = x[5] ^ x[21]
eor r4, r4, r3    @ x[21] = x[21] ^ x[11]
eor r9, r9, r6    @ x[5] = x[5] ^ x[7]
eor r6, r6, r7    @ x[7] = x[7] ^ x[23]
eor r3, r3, r8    @ x[11] = x[11] ^ x[22]
vmov s5, r9    @ spill x[5] from r9
vmov r9, s27    @ load x[27] from s27
eor r8, r8, r9    @ x[22] = x[22] ^ x[27]
eor r8, r8, r7    @ x[22] = x[22] ^ x[23]
eor r7, r7, r10    @ x[23] = x[23] ^ x[20]
eor r10, r10, r9    @ x[20] = x[20] ^ x[27]
vmov s22, r8    @ spill x[22] from r8
vmov r8, s15    @ load x[15] from s15
eor r9, r9, r8    @ x[27] = x[27] ^ x[15]
eor r6, r6, r1    @ x[7] = x[7] ^ x[1]
eor r8, r8, r1    @ x[15] = x[15] ^ x[1]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s19    @ load x[19] from s19
eor r1, r1, r6    @ x[1] = x[1] ^ x[19]
vmov s20, r10    @ spill x[20] from r10
vmov r10, s28    @ load x[28] from s28
eor r6, r6, r10    @ x[19] = x[19] ^ x[28]
eor r1, r1, r2    @ x[1] = x[1] ^ x[25]
eor r10, r10, r2    @ x[28] = x[28] ^ x[25]
vmov s1, r1    @ spill x[1] from r1
vmov r1, s4    @ load x[4] from s4
vmov s11, r3    @ spill x[11] from r3
vmov r3, s10    @ load x[10] from s10
eor r1, r1, r3    @ x[4] = x[4] ^ x[10]
eor r5, r5, r11    @ x[31] = x[31] ^ x[16]
vmov s31, r5    @ spill x[31] from r5
vmov r5, s2    @ load x[2] from s2
eor r2, r2, r5    @ x[25] = x[25] ^ x[2]
vmov s25, r2    @ spill x[25] from r2
vmov r2, s29    @ load x[29] from s29
eor r5, r5, r2    @ x[2] = x[2] ^ x[29]
eor r5, r5, r8    @ x[2] = x[2] ^ x[15]
vmov s10, r3    @ spill x[10] from r3
vmov r3, s26    @ load x[26] from s26
eor r6, r6, r3    @ x[19] = x[19] ^ x[26]
eor r9, r9, r0    @ x[27] = x[27] ^ x[14]
vmov s27, r9    @ spill x[27] from r9
vmov r9, s13    @ load x[13] from s13
eor r6, r6, r9    @ x[19] = x[19] ^ x[13]
eor r11, r11, r10    @ x[16] = x[16] ^ x[28]
eor r10, r10, r1    @ x[28] = x[28] ^ x[4]
vmov s16, r11    @ spill x[16] from r11
vmov r11, s6    @ load x[6] from s6
eor r1, r1, r11    @ x[4] = x[4] ^ x[6]
eor r10, r10, r0    @ x[28] = x[28] ^ x[14]
eor r2, r2, r4    @ x[29] = x[29] ^ x[21]
vmov s29, r2    @ spill x[29] from r2
vmov r2, s24    @ load x[24] from s24
eor r6, r6, r2    @ x[19] = x[19] ^ x[24]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s8    @ load x[8] from s8
eor r8, r8, r6    @ x[15] = x[15] ^ x[8]
eor r7, r7, r2    @ x[23] = x[23] ^ x[24]
eor r2, r2, r1    @ x[24] = x[24] ^ x[4]
vmov s28, r10    @ spill x[28] from r10
vmov r10, s9    @ load x[9] from s9
eor r2, r2, r10    @ x[24] = x[24] ^ x[9]
vmov s24, r2    @ spill x[24] from r2
vmov r2, s11    @ load x[11] from s11
eor r10, r10, r2    @ x[9] = x[9] ^ x[11]
vmov s9, r10    @ spill x[9] from r10
vmov r10, s20    @ load x[20] from s20
eor r10, r10, r3    @ x[20] = x[20] ^ x[26]
vmov s23, r7    @ spill x[23] from r7
vmov r7, s30    @ load x[30] from s30
eor r11, r11, r7    @ x[6] = x[6] ^ x[30]
eor r3, r3, r7    @ x[26] = x[26] ^ x[30]
vmov s15, r8    @ spill x[15] from r8
vmov r8, s7    @ load x[7] from s7
eor r7, r7, r8    @ x[30] = x[30] ^ x[7]
eor r8, r8, r5    @ x[7] = x[7] ^ x[2]
eor r5, r5, r10    @ x[2] = x[2] ^ x[20]
eor r6, r6, r3    @ x[8] = x[8] ^ x[26]
eor r2, r2, r3    @ x[11] = x[11] ^ x[26]
eor r5, r5, r6    @ x[2] = x[2] ^ x[8]
eor r3, r3, r11    @ x[26] = x[26] ^ x[6]
vmov s2, r5    @ spill x[2] from r5
vmov r5, s22    @ load x[22] from s22
eor r8, r8, r5    @ x[7] = x[7] ^ x[22]
eor r5, r5, r2    @ x[22] = x[22] ^ x[11]
eor r6, r6, r5    @ x[8] = x[8] ^ x[22]
eor r1, r1, r4    @ x[4] = x[4] ^ x[21]
eor r2, r2, r9    @ x[11] = x[11] ^ x[13]
vmov s11, r2    @ spill x[11] from r2
vmov r2, s12    @ load x[12] from s12
eor r6, r6, r2    @ x[8] = x[8] ^ x[12]
eor r0, r0, r2    @ x[14] = x[14] ^ x[12]
vmov s8, r6    @ spill x[8] from r6
vmov r6, s18    @ load x[18] from s18
eor r2, r2, r6    @ x[12] = x[12] ^ x[18]
vmov s7, r8    @ spill x[7] from r8
vmov r8, s17    @ load x[17] from s17
eor r9, r9, r8    @ x[13] = x[13] ^ x[17]
vmov s30, r7    @ spill x[30] from r7
vmov r7, s15    @ load x[15] from s15
vmov s13, r9    @ spill x[13] from r9
vmov r9, s10    @ load x[10] from s10
eor r7, r7, r9    @ x[15] = x[15] ^ x[10]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s16    @ load x[16] from s16
eor r7, r7, r4    @ x[16] = x[16] ^ x[21]
vmov s16, r7    @ spill x[16] from r7
vmov r7, s31    @ load x[31] from s31
eor r4, r4, r7    @ x[21] = x[21] ^ x[31]
vmov s12, r2    @ spill x[12] from r2
vmov r2, s25    @ load x[25] from s25
eor r5, r5, r2    @ x[22] = x[22] ^ x[25]
eor r6, r6, r7    @ x[18] = x[18] ^ x[31]
eor r6, r6, r5    @ x[18] = x[18] ^ x[22]
eor r5, r5, r9    @ x[22] = x[22] ^ x[10]
eor r9, r9, r8    @ x[10] = x[10] ^ x[17]
vmov s22, r5    @ spill x[22] from r5
vmov r5, s1    @ load x[1] from s1
vmov s31, r7    @ spill x[31] from r7
vmov r7, s3    @ load x[3] from s3
eor r5, r5, r7    @ x[1] = x[1] ^ x[3]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s27    @ load x[27] from s27
eor r6, r6, r3    @ x[27] = x[27] ^ x[26]
vmov s27, r6    @ spill x[27] from r6
vmov r6, s23    @ load x[23] from s23
eor r3, r3, r6    @ x[26] = x[26] ^ x[23]
eor r8, r8, r11    @ x[17] = x[17] ^ x[6]
eor r11, r11, r1    @ x[6] = x[6] ^ x[4]
eor r0, r0, r7    @ x[14] = x[14] ^ x[3]
vmov s17, r8    @ spill x[17] from r8
vmov r8, s28    @ load x[28] from s28
vmov s4, r1    @ spill x[4] from r1
vmov r1, s0    @ load x[0] from s0
eor r8, r8, r1    @ x[28] = x[28] ^ x[0]
eor r4, r4, r10    @ x[21] = x[21] ^ x[20]
vmov s21, r4    @ spill x[21] from r4
vmov r4, s13    @ load x[13] from s13
vmov s20, r10    @ spill x[20] from r10
vmov r10, s12    @ load x[12] from s12
eor r4, r4, r10    @ x[13] = x[13] ^ x[12]
vmov s13, r4    @ spill x[13] from r4
vmov r4, s30    @ load x[30] from s30
eor r1, r1, r4    @ x[0] = x[0] ^ x[30]
vmov s0, r1    @ spill x[0] from r1
vmov r1, s29    @ load x[29] from s29
eor r9, r9, r1    @ x[10] = x[10] ^ x[29]
eor r2, r2, r0    @ x[25] = x[25] ^ x[14]
vmov s29, r1    @ spill x[29] from r1
vmov r1, s18    @ load x[18] from s18
eor r1, r1, r5    @ x[18] = x[18] ^ x[1]
eor r6, r6, r5    @ x[23] = x[23] ^ x[1]
eor r2, r2, r8    @ x[25] = x[25] ^ x[28]
eor r10, r10, r11    @ x[12] = x[12] ^ x[6]
vmov s23, r6    @ spill x[23] from r6
vmov r6, s7    @ load x[7] from s7
vmov s12, r10    @ spill x[12] from r10
vmov r10, s16    @ load x[16] from s16
eor r6, r6, r10    @ x[7] = x[7] ^ x[16]
vmov s16, r10    @ spill x[16] from r10
vmov r10, s2    @ load x[2] from s2
eor r0, r0, r10    @ x[14] = x[14] ^ x[2]
eor r11, r11, r10    @ x[6] = x[6] ^ x[2]
vmov s14, r0    @ spill x[14] from r0
vmov r0, s31    @ load x[31] from s31
eor r10, r10, r0    @ x[2] = x[2] ^ x[31]
vmov s18, r1    @ spill x[18] from r1
vmov r1, s9    @ load x[9] from s9
eor r1, r1, r9    @ x[9] = x[9] ^ x[10]
eor r0, r0, r1    @ x[31] = x[31] ^ x[9]
eor r1, r1, r5    @ x[9] = x[9] ^ x[1]
eor r11, r11, r3    @ x[6] = x[6] ^ x[26]
vmov s6, r11    @ spill x[6] from r11
vmov r11, s5    @ load x[5] from s5
eor r0, r0, r11    @ x[31] = x[31] ^ x[5]
eor r11, r11, r7    @ x[5] = x[5] ^ x[3]
eor r3, r3, r7    @ x[26] = x[26] ^ x[3]
eor r5, r5, r7    @ x[1] = x[1] ^ x[3]
vmov s5, r11    @ spill x[5] from r11
vmov r11, s24    @ load x[24] from s24
eor r7, r7, r11    @ x[3] = x[3] ^ x[24]
eor r0, r0, r11    @ x[31] = x[31] ^ x[24]
vmov s31, r0    @ spill x[31] from r0
vmov r0, s8    @ load x[8] from s8
eor r11, r11, r0    @ x[24] = x[24] ^ x[8]
eor r4, r4, r3    @ x[30] = x[30] ^ x[26]
eor r3, r3, r0    @ x[26] = x[26] ^ x[8]
vmov s30, r4    @ spill x[30] from r4
vmov r4, s4    @ load x[4] from s4
vmov s26, r3    @ spill x[26] from r3
vmov r3, s29    @ load x[29] from s29
eor r4, r4, r3    @ x[4] = x[4] ^ x[29]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s15    @ load x[15] from s15
eor r4, r4, r10    @ x[4] = x[4] ^ x[15]
eor r1, r1, r10    @ x[9] = x[9] ^ x[15]
eor r5, r5, r3    @ x[1] = x[1] ^ x[29]
eor r10, r10, r2    @ x[15] = x[15] ^ x[25]
vmov s15, r10    @ spill x[15] from r10
vmov r10, s17    @ load x[17] from s17
vmov s9, r1    @ spill x[9] from r1
vmov r1, s20    @ load x[20] from s20
eor r10, r10, r1    @ x[17] = x[17] ^ x[20]
vmov s20, r1    @ spill x[20] from r1
vmov r1, s27    @ load x[27] from s27
eor r0, r0, r1    @ x[8] = x[8] ^ x[27]
eor r3, r3, r6    @ x[29] = x[29] ^ x[7]
eor r0, r0, r2    @ x[8] = x[8] ^ x[25]
vmov s25, r2    @ spill x[25] from r2
vmov r2, s22    @ load x[22] from s22
eor r11, r11, r2    @ x[24] = x[24] ^ x[22]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s19    @ load x[19] from s19
eor r0, r0, r6    @ x[8] = x[8] ^ x[19]
eor r10, r10, r1    @ x[17] = x[17] ^ x[27]
eor r4, r4, r8    @ x[4] = x[4] ^ x[28]
eor r2, r2, r9    @ x[22] = x[22] ^ x[10]
vmov s10, r9    @ spill x[10] from r9
vmov r9, s11    @ load x[11] from s11
eor r2, r2, r9    @ x[22] = x[22] ^ x[11]
eor r1, r1, r7    @ x[27] = x[27] ^ x[3]
vmov s8, r0    @ spill x[8] from r0
vmov r0, s2    @ load x[2] from s2
eor r3, r3, r0    @ x[29] = x[29] ^ x[2]
vmov s4, r4    @ spill x[4] from r4
vmov r4, s0    @ load x[0] from s0
eor r5, r5, r4    @ x[1] = x[1] ^ x[0]
eor r7, r7, r9    @ x[3] = x[3] ^ x[11]
vmov s11, r9    @ spill x[11] from r9
vmov r9, s13    @ load x[13] from s13
eor r4, r4, r9    @ x[0] = x[0] ^ x[13]
vmov s17, r10    @ spill x[17] from r10
vmov r10, s26    @ load x[26] from s26
eor r10, r10, r5    @ x[26] = x[26] ^ x[1]
eor r8, r8, r7    @ x[28] = x[28] ^ x[3]
eor r0, r0, r1    @ x[2] = x[2] ^ x[27]
eor r6, r6, r2    @ x[19] = x[19] ^ x[22]
vmov s1, r5    @ spill x[1] from r5
vmov r5, s31    @ load x[31] from s31
vmov s3, r7    @ spill x[3] from r7
vmov r7, s30    @ load x[30] from s30
eor r5, r5, r7    @ x[31] = x[31] ^ x[30]
eor r0, r0, r4    @ x[2] = x[2] ^ x[0]
vmov s28, r8    @ spill x[28] from r8
vmov r8, s9    @ load x[9] from s9
eor r3, r3, r8    @ x[29] = x[29] ^ x[9]
vmov s29, r3    @ spill x[29] from r3
vmov r3, s21    @ load x[21] from s21
eor r3, r3, r11    @ x[21] = x[21] ^ x[24]
vmov s24, r11    @ spill x[24] from r11
vmov r11, s18    @ load x[18] from s18
eor r9, r9, r11    @ x[13] = x[13] ^ x[18]
vmov s27, r1    @ spill x[27] from r1
vmov r1, s6    @ load x[6] from s6
eor r1, r1, r0    @ x[6] = x[6] ^ x[2]
vmov s22, r2    @ spill x[22] from r2
vmov r2, s17    @ load x[17] from s17
eor r2, r2, r6    @ x[17] = x[17] ^ x[19]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s14    @ load x[14] from s14
eor r6, r6, r0    @ x[14] = x[14] ^ x[2]
vstr.32  s15, [r14, #0]    @ y[0] = x[15] from s15
vstr.32  s28, [r14, #4]    @ y[1] = x[28] from s28
vstr.32  s25, [r14, #8]    @ y[2] = x[25] from s25
vstr.32  s11, [r14, #12]    @ y[3] = x[11] from s11
vstr.32  s4, [r14, #16]    @ y[4] = x[4] from s4
vstr.32  s27, [r14, #20]    @ y[5] = x[27] from s27
vstr.32  s3, [r14, #24]    @ y[6] = x[3] from s3
str  r0, [r14, #28]    @ y[7] = x[2] from r0
str  r1, [r14, #32]    @ y[8] = x[6] from r1
str  r3, [r14, #36]    @ y[9] = x[21] from r3
str  r4, [r14, #40]    @ y[10] = x[0] from r4
vstr.32  s12, [r14, #44]    @ y[11] = x[12] from s12
vstr.32  s24, [r14, #48]    @ y[12] = x[24] from s24
vstr.32  s8, [r14, #52]    @ y[13] = x[8] from s8
str  r10, [r14, #56]    @ y[14] = x[26] from r10
str  r6, [r14, #60]    @ y[15] = x[14] from r6
vstr.32  s23, [r14, #64]    @ y[16] = x[23] from s23
vstr.32  s1, [r14, #68]    @ y[17] = x[1] from s1
vstr.32  s20, [r14, #72]    @ y[18] = x[20] from s20
vstr.32  s5, [r14, #76]    @ y[19] = x[5] from s5
str  r11, [r14, #80]    @ y[20] = x[18] from r11
str  r5, [r14, #84]    @ y[21] = x[31] from r5
vstr.32  s19, [r14, #88]    @ y[22] = x[19] from s19
str  r7, [r14, #92]    @ y[23] = x[30] from r7
vstr.32  s29, [r14, #96]    @ y[24] = x[29] from s29
str  r8, [r14, #100]    @ y[25] = x[9] from r8
vstr.32  s16, [r14, #104]    @ y[26] = x[16] from s16
vstr.32  s10, [r14, #108]    @ y[27] = x[10] from s10
vstr.32  s22, [r14, #112]    @ y[28] = x[22] from s22
str  r9, [r14, #116]    @ y[29] = x[13] from r9
vstr.32  s7, [r14, #120]    @ y[30] = x[7] from s7
str  r2, [r14, #124]    @ y[31] = x[17] from r2
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v17, .-gft_mul_v17
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v18
.type gft_mul_v18, %function
.align 2
gft_mul_v18:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s8    @ load x[8] from s8
vmov r1, s24    @ load x[24] from s24
eor r0, r0, r1    @ x[8] = x[8] ^ x[24]
vmov r2, s11    @ load x[11] from s11
eor r2, r2, r0    @ x[11] = x[11] ^ x[8]
vmov r3, s10    @ load x[10] from s10
eor r3, r3, r2    @ x[10] = x[10] ^ x[11]
vmov r4, s26    @ load x[26] from s26
vmov r5, s0    @ load x[0] from s0
eor r4, r4, r5    @ x[26] = x[26] ^ x[0]
vmov r6, s14    @ load x[14] from s14
vmov r7, s1    @ load x[1] from s1
eor r6, r6, r7    @ x[14] = x[14] ^ x[1]
vmov r8, s30    @ load x[30] from s30
eor r0, r0, r8    @ x[8] = x[8] ^ x[30]
vmov r9, s13    @ load x[13] from s13
vmov r10, s28    @ load x[28] from s28
eor r9, r9, r10    @ x[13] = x[13] ^ x[28]
vmov r11, s23    @ load x[23] from s23
eor r0, r0, r11    @ x[8] = x[8] ^ x[23]
eor r11, r11, r3    @ x[23] = x[23] ^ x[10]
vmov s8, r0    @ spill x[8] from r0
vmov r0, s18    @ load x[18] from s18
vmov s24, r1    @ spill x[24] from r1
vmov r1, s22    @ load x[22] from s22
eor r0, r0, r1    @ x[18] = x[18] ^ x[22]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s19    @ load x[19] from s19
eor r6, r6, r10    @ x[19] = x[19] ^ x[28]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s12    @ load x[12] from s12
eor r0, r0, r6    @ x[18] = x[18] ^ x[12]
vmov s18, r0    @ spill x[18] from r0
vmov r0, s9    @ load x[9] from s9
eor r0, r0, r1    @ x[9] = x[9] ^ x[22]
eor r9, r9, r3    @ x[13] = x[13] ^ x[10]
vmov s13, r9    @ spill x[13] from r9
vmov r9, s5    @ load x[5] from s5
vmov s9, r0    @ spill x[9] from r0
vmov r0, s6    @ load x[6] from s6
eor r9, r9, r0    @ x[5] = x[5] ^ x[6]
vmov s10, r3    @ spill x[10] from r3
vmov r3, s21    @ load x[21] from s21
eor r0, r0, r3    @ x[6] = x[6] ^ x[21]
vmov s6, r0    @ spill x[6] from r0
vmov r0, s20    @ load x[20] from s20
eor r0, r0, r1    @ x[20] = x[20] ^ x[22]
vmov s22, r1    @ spill x[22] from r1
vmov r1, s29    @ load x[29] from s29
eor r2, r2, r1    @ x[11] = x[11] ^ x[29]
eor r1, r1, r8    @ x[29] = x[29] ^ x[30]
vmov s11, r2    @ spill x[11] from r2
vmov r2, s7    @ load x[7] from s7
eor r0, r0, r2    @ x[20] = x[20] ^ x[7]
eor r9, r9, r11    @ x[5] = x[5] ^ x[23]
eor r8, r8, r10    @ x[30] = x[30] ^ x[28]
vmov s5, r9    @ spill x[5] from r9
vmov r9, s3    @ load x[3] from s3
eor r5, r5, r9    @ x[0] = x[0] ^ x[3]
eor r3, r3, r2    @ x[21] = x[21] ^ x[7]
vmov s20, r0    @ spill x[20] from r0
vmov r0, s2    @ load x[2] from s2
eor r0, r0, r4    @ x[2] = x[2] ^ x[26]
eor r10, r10, r7    @ x[28] = x[28] ^ x[1]
eor r4, r4, r7    @ x[26] = x[26] ^ x[1]
eor r9, r9, r6    @ x[3] = x[3] ^ x[12]
eor r6, r6, r1    @ x[12] = x[12] ^ x[29]
eor r4, r4, r11    @ x[26] = x[26] ^ x[23]
eor r4, r4, r6    @ x[26] = x[26] ^ x[12]
vmov s26, r4    @ spill x[26] from r4
vmov r4, s31    @ load x[31] from s31
eor r6, r6, r4    @ x[12] = x[12] ^ x[31]
eor r7, r7, r2    @ x[1] = x[1] ^ x[7]
vmov s30, r8    @ spill x[30] from r8
vmov r8, s6    @ load x[6] from s6
eor r4, r4, r8    @ x[31] = x[31] ^ x[6]
vmov s12, r6    @ spill x[12] from r6
vmov r6, s11    @ load x[11] from s11
vmov s21, r3    @ spill x[21] from r3
vmov r3, s22    @ load x[22] from s22
eor r6, r6, r3    @ x[11] = x[11] ^ x[22]
eor r11, r11, r2    @ x[23] = x[23] ^ x[7]
vmov s23, r11    @ spill x[23] from r11
vmov r11, s10    @ load x[10] from s10
eor r11, r11, r3    @ x[10] = x[10] ^ x[22]
vmov s11, r6    @ spill x[11] from r6
vmov r6, s18    @ load x[18] from s18
eor r6, r6, r8    @ x[18] = x[18] ^ x[6]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s19    @ load x[19] from s19
eor r6, r6, r1    @ x[19] = x[19] ^ x[29]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s9    @ load x[9] from s9
eor r2, r2, r6    @ x[7] = x[7] ^ x[9]
eor r3, r3, r5    @ x[22] = x[22] ^ x[0]
eor r8, r8, r11    @ x[6] = x[6] ^ x[10]
eor r5, r5, r8    @ x[0] = x[0] ^ x[6]
eor r1, r1, r8    @ x[29] = x[29] ^ x[6]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s4    @ load x[4] from s4
eor r8, r8, r11    @ x[6] = x[6] ^ x[4]
eor r6, r6, r11    @ x[9] = x[9] ^ x[4]
eor r10, r10, r4    @ x[28] = x[28] ^ x[31]
eor r10, r10, r1    @ x[28] = x[28] ^ x[29]
eor r7, r7, r3    @ x[1] = x[1] ^ x[22]
vmov s29, r1    @ spill x[29] from r1
vmov r1, s27    @ load x[27] from s27
vmov s28, r10    @ spill x[28] from r10
vmov r10, s13    @ load x[13] from s13
eor r1, r1, r10    @ x[27] = x[27] ^ x[13]
eor r11, r11, r0    @ x[4] = x[4] ^ x[2]
eor r3, r3, r1    @ x[22] = x[22] ^ x[27]
eor r7, r7, r9    @ x[1] = x[1] ^ x[3]
vmov s27, r1    @ spill x[27] from r1
vmov r1, s14    @ load x[14] from s14
eor r5, r5, r1    @ x[0] = x[0] ^ x[14]
vmov s0, r5    @ spill x[0] from r5
vmov r5, s24    @ load x[24] from s24
eor r0, r0, r5    @ x[2] = x[2] ^ x[24]
eor r1, r1, r5    @ x[14] = x[14] ^ x[24]
eor r9, r9, r0    @ x[3] = x[3] ^ x[2]
eor r0, r0, r2    @ x[2] = x[2] ^ x[7]
eor r2, r2, r4    @ x[7] = x[7] ^ x[31]
eor r4, r4, r8    @ x[31] = x[31] ^ x[6]
vmov s7, r2    @ spill x[7] from r2
vmov r2, s10    @ load x[10] from s10
vmov s13, r10    @ spill x[13] from r10
vmov r10, s11    @ load x[11] from s11
eor r2, r2, r10    @ x[10] = x[10] ^ x[11]
vmov s4, r11    @ spill x[4] from r11
vmov r11, s21    @ load x[21] from s21
eor r3, r3, r11    @ x[22] = x[22] ^ x[21]
vmov s14, r1    @ spill x[14] from r1
vmov r1, s12    @ load x[12] from s12
eor r1, r1, r5    @ x[12] = x[12] ^ x[24]
vmov s21, r11    @ spill x[21] from r11
vmov r11, s18    @ load x[18] from s18
eor r11, r11, r1    @ x[18] = x[18] ^ x[12]
vmov s18, r11    @ spill x[18] from r11
vmov r11, s30    @ load x[30] from s30
eor r4, r4, r11    @ x[31] = x[31] ^ x[30]
eor r6, r6, r4    @ x[9] = x[9] ^ x[31]
eor r8, r8, r9    @ x[6] = x[6] ^ x[3]
eor r0, r0, r2    @ x[2] = x[2] ^ x[10]
vmov s3, r9    @ spill x[3] from r9
vmov r9, s17    @ load x[17] from s17
eor r10, r10, r9    @ x[11] = x[11] ^ x[17]
eor r0, r0, r3    @ x[2] = x[2] ^ x[22]
eor r10, r10, r11    @ x[11] = x[11] ^ x[30]
vmov s2, r0    @ spill x[2] from r0
vmov r0, s23    @ load x[23] from s23
eor r4, r4, r0    @ x[31] = x[31] ^ x[23]
eor r8, r8, r1    @ x[6] = x[6] ^ x[12]
vmov s17, r9    @ spill x[17] from r9
vmov r9, s15    @ load x[15] from s15
eor r11, r11, r9    @ x[30] = x[30] ^ x[15]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s20    @ load x[20] from s20
eor r4, r4, r11    @ x[31] = x[31] ^ x[20]
eor r8, r8, r6    @ x[6] = x[6] ^ x[9]
eor r6, r6, r7    @ x[9] = x[9] ^ x[1]
vmov s6, r8    @ spill x[6] from r8
vmov r8, s16    @ load x[16] from s16
eor r4, r4, r8    @ x[31] = x[31] ^ x[16]
eor r8, r8, r7    @ x[16] = x[16] ^ x[1]
eor r10, r10, r5    @ x[11] = x[11] ^ x[24]
vmov s11, r10    @ spill x[11] from r10
vmov r10, s8    @ load x[8] from s8
eor r6, r6, r10    @ x[9] = x[9] ^ x[8]
vmov s31, r4    @ spill x[31] from r4
vmov r4, s21    @ load x[21] from s21
eor r7, r7, r4    @ x[1] = x[1] ^ x[21]
vmov s9, r6    @ spill x[9] from r6
vmov r6, s19    @ load x[19] from s19
eor r4, r4, r6    @ x[21] = x[21] ^ x[19]
eor r6, r6, r3    @ x[19] = x[19] ^ x[22]
vmov s21, r4    @ spill x[21] from r4
vmov r4, s25    @ load x[25] from s25
eor r7, r7, r4    @ x[1] = x[1] ^ x[25]
eor r9, r9, r0    @ x[15] = x[15] ^ x[23]
eor r2, r2, r9    @ x[10] = x[10] ^ x[15]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s14    @ load x[14] from s14
eor r0, r0, r6    @ x[23] = x[23] ^ x[14]
vmov s16, r8    @ spill x[16] from r8
vmov r8, s4    @ load x[4] from s4
eor r6, r6, r8    @ x[14] = x[14] ^ x[4]
eor r6, r6, r1    @ x[14] = x[14] ^ x[12]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s13    @ load x[13] from s13
eor r1, r1, r6    @ x[12] = x[12] ^ x[13]
eor r8, r8, r10    @ x[4] = x[4] ^ x[8]
eor r6, r6, r11    @ x[13] = x[13] ^ x[20]
eor r11, r11, r5    @ x[20] = x[20] ^ x[24]
vmov s13, r6    @ spill x[13] from r6
vmov r6, s3    @ load x[3] from s3
eor r5, r5, r6    @ x[24] = x[24] ^ x[3]
eor r10, r10, r11    @ x[8] = x[8] ^ x[20]
eor r9, r9, r3    @ x[15] = x[15] ^ x[22]
eor r4, r4, r3    @ x[25] = x[25] ^ x[22]
eor r11, r11, r3    @ x[20] = x[20] ^ x[22]
vmov s8, r10    @ spill x[8] from r10
vmov r10, s27    @ load x[27] from s27
eor r3, r3, r10    @ x[22] = x[22] ^ x[27]
vmov s20, r11    @ spill x[20] from r11
vmov r11, s28    @ load x[28] from s28
eor r10, r10, r11    @ x[27] = x[27] ^ x[28]
vmov s27, r10    @ spill x[27] from r10
vmov r10, s5    @ load x[5] from s5
eor r8, r8, r10    @ x[4] = x[4] ^ x[5]
vmov s4, r8    @ spill x[4] from r8
vmov r8, s26    @ load x[26] from s26
eor r11, r11, r8    @ x[28] = x[28] ^ x[26]
eor r1, r1, r3    @ x[12] = x[12] ^ x[22]
eor r3, r3, r6    @ x[22] = x[22] ^ x[3]
vmov s22, r3    @ spill x[22] from r3
vmov r3, s29    @ load x[29] from s29
eor r10, r10, r3    @ x[5] = x[5] ^ x[29]
eor r5, r5, r3    @ x[24] = x[24] ^ x[29]
vmov s3, r6    @ spill x[3] from r6
vmov r6, s0    @ load x[0] from s0
eor r2, r2, r6    @ x[10] = x[10] ^ x[0]
eor r7, r7, r6    @ x[1] = x[1] ^ x[0]
eor r9, r9, r10    @ x[15] = x[15] ^ x[5]
vmov s15, r9    @ spill x[15] from r9
vmov r9, s17    @ load x[17] from s17
eor r10, r10, r9    @ x[5] = x[5] ^ x[17]
eor r6, r6, r4    @ x[0] = x[0] ^ x[25]
eor r4, r4, r9    @ x[25] = x[25] ^ x[17]
eor r8, r8, r0    @ x[26] = x[26] ^ x[23]
eor r7, r7, r8    @ x[1] = x[1] ^ x[26]
vmov s0, r6    @ spill x[0] from r6
vmov r6, s16    @ load x[16] from s16
eor r9, r9, r6    @ x[17] = x[17] ^ x[16]
vmov s17, r9    @ spill x[17] from r9
vmov r9, s30    @ load x[30] from s30
eor r8, r8, r9    @ x[26] = x[26] ^ x[30]
vmov s26, r8    @ spill x[26] from r8
vmov r8, s18    @ load x[18] from s18
eor r0, r0, r8    @ x[23] = x[23] ^ x[18]
vmov s1, r7    @ spill x[1] from r7
vmov r7, s6    @ load x[6] from s6
eor r3, r3, r7    @ x[29] = x[29] ^ x[6]
vmov s29, r3    @ spill x[29] from r3
vmov r3, s19    @ load x[19] from s19
eor r8, r8, r3    @ x[18] = x[18] ^ x[19]
eor r5, r5, r0    @ x[24] = x[24] ^ x[23]
vmov s24, r5    @ spill x[24] from r5
vmov r5, s20    @ load x[20] from s20
eor r5, r5, r8    @ x[20] = x[20] ^ x[18]
vmov s18, r8    @ spill x[18] from r8
vmov r8, s8    @ load x[8] from s8
vmov s5, r10    @ spill x[5] from r10
vmov r10, s21    @ load x[21] from s21
eor r8, r8, r10    @ x[8] = x[8] ^ x[21]
eor r7, r7, r11    @ x[6] = x[6] ^ x[28]
eor r3, r3, r1    @ x[19] = x[19] ^ x[12]
eor r1, r1, r9    @ x[12] = x[12] ^ x[30]
vmov s6, r7    @ spill x[6] from r7
vmov r7, s7    @ load x[7] from s7
eor r7, r7, r5    @ x[7] = x[7] ^ x[20]
vmov s8, r8    @ spill x[8] from r8
vmov r8, s3    @ load x[3] from s3
vmov s7, r7    @ spill x[7] from r7
vmov r7, s13    @ load x[13] from s13
eor r8, r8, r7    @ x[3] = x[3] ^ x[13]
eor r0, r0, r2    @ x[23] = x[23] ^ x[10]
vmov s23, r0    @ spill x[23] from r0
vmov r0, s22    @ load x[22] from s22
eor r0, r0, r4    @ x[22] = x[22] ^ x[25]
eor r2, r2, r1    @ x[10] = x[10] ^ x[12]
eor r6, r6, r7    @ x[16] = x[16] ^ x[13]
eor r1, r1, r4    @ x[12] = x[12] ^ x[25]
eor r4, r4, r7    @ x[25] = x[25] ^ x[13]
eor r7, r7, r5    @ x[13] = x[13] ^ x[20]
vmov s10, r2    @ spill x[10] from r2
vmov r2, s5    @ load x[5] from s5
eor r5, r5, r2    @ x[20] = x[20] ^ x[5]
vmov s12, r1    @ spill x[12] from r1
vmov r1, s9    @ load x[9] from s9
eor r9, r9, r1    @ x[30] = x[30] ^ x[9]
vmov s25, r4    @ spill x[25] from r4
vmov r4, s18    @ load x[18] from s18
vmov s20, r5    @ spill x[20] from r5
vmov r5, s1    @ load x[1] from s1
eor r4, r4, r5    @ x[18] = x[18] ^ x[1]
vmov s5, r2    @ spill x[5] from r2
vmov r2, s14    @ load x[14] from s14
vmov s1, r5    @ spill x[1] from r5
vmov r5, s15    @ load x[15] from s15
eor r2, r2, r5    @ x[14] = x[14] ^ x[15]
vmov s14, r2    @ spill x[14] from r2
vmov r2, s17    @ load x[17] from s17
eor r1, r1, r2    @ x[9] = x[9] ^ x[17]
vmov s9, r1    @ spill x[9] from r1
vmov r1, s27    @ load x[27] from s27
eor r2, r2, r1    @ x[17] = x[17] ^ x[27]
vmov s15, r5    @ spill x[15] from r5
vmov r5, s26    @ load x[26] from s26
eor r2, r2, r5    @ x[17] = x[17] ^ x[26]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s4    @ load x[4] from s4
eor r8, r8, r2    @ x[3] = x[3] ^ x[4]
vmov s3, r8    @ spill x[3] from r8
vmov r8, s31    @ load x[31] from s31
eor r10, r10, r8    @ x[21] = x[21] ^ x[31]
vmov s21, r10    @ spill x[21] from r10
vmov r10, s2    @ load x[2] from s2
eor r7, r7, r10    @ x[13] = x[13] ^ x[2]
eor r1, r1, r4    @ x[27] = x[27] ^ x[18]
eor r7, r7, r1    @ x[13] = x[13] ^ x[27]
eor r0, r0, r6    @ x[22] = x[22] ^ x[16]
vmov s16, r6    @ spill x[16] from r6
vmov r6, s24    @ load x[24] from s24
eor r6, r6, r9    @ x[24] = x[24] ^ x[30]
eor r11, r11, r3    @ x[28] = x[28] ^ x[19]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s23    @ load x[23] from s23
eor r9, r9, r11    @ x[23] = x[23] ^ x[28]
eor r10, r10, r0    @ x[2] = x[2] ^ x[22]
eor r2, r2, r7    @ x[4] = x[4] ^ x[13]
str  r2, [r14, #0]    @ y[0] = x[4] from r2
vstr.32  s29, [r14, #4]    @ y[1] = x[29] from s29
vstr.32  s1, [r14, #8]    @ y[2] = x[1] from s1
str  r7, [r14, #12]    @ y[3] = x[13] from r7
vstr.32  s15, [r14, #16]    @ y[4] = x[15] from s15
vstr.32  s14, [r14, #20]    @ y[5] = x[14] from s14
str  r5, [r14, #24]    @ y[6] = x[26] from r5
vstr.32  s0, [r14, #28]    @ y[7] = x[0] from s0
vstr.32  s10, [r14, #32]    @ y[8] = x[10] from s10
vstr.32  s30, [r14, #36]    @ y[9] = x[30] from s30
vstr.32  s3, [r14, #40]    @ y[10] = x[3] from s3
vstr.32  s6, [r14, #44]    @ y[11] = x[6] from s6
str  r10, [r14, #48]    @ y[12] = x[2] from r10
vstr.32  s16, [r14, #52]    @ y[13] = x[16] from s16
str  r9, [r14, #56]    @ y[14] = x[23] from r9
str  r6, [r14, #60]    @ y[15] = x[24] from r6
vstr.32  s25, [r14, #64]    @ y[16] = x[25] from s25
str  r0, [r14, #68]    @ y[17] = x[22] from r0
str  r11, [r14, #72]    @ y[18] = x[28] from r11
str  r1, [r14, #76]    @ y[19] = x[27] from r1
vstr.32  s21, [r14, #80]    @ y[20] = x[21] from s21
vstr.32  s20, [r14, #84]    @ y[21] = x[20] from s20
str  r8, [r14, #88]    @ y[22] = x[31] from r8
vstr.32  s5, [r14, #92]    @ y[23] = x[5] from s5
vstr.32  s11, [r14, #96]    @ y[24] = x[11] from s11
vstr.32  s9, [r14, #100]    @ y[25] = x[9] from s9
vstr.32  s8, [r14, #104]    @ y[26] = x[8] from s8
vstr.32  s7, [r14, #108]    @ y[27] = x[7] from s7
vstr.32  s12, [r14, #112]    @ y[28] = x[12] from s12
vstr.32  s17, [r14, #116]    @ y[29] = x[17] from s17
str  r3, [r14, #120]    @ y[30] = x[19] from r3
str  r4, [r14, #124]    @ y[31] = x[18] from r4
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v18, .-gft_mul_v18
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v19
.type gft_mul_v19, %function
.align 2
gft_mul_v19:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s15    @ load x[15] from s15
vmov r1, s13    @ load x[13] from s13
eor r0, r0, r1    @ x[15] = x[15] ^ x[13]
vmov r2, s5    @ load x[5] from s5
vmov r3, s28    @ load x[28] from s28
eor r2, r2, r3    @ x[5] = x[5] ^ x[28]
vmov r4, s6    @ load x[6] from s6
eor r4, r4, r2    @ x[6] = x[6] ^ x[5]
eor r4, r4, r0    @ x[6] = x[6] ^ x[15]
vmov r5, s30    @ load x[30] from s30
vmov r6, s26    @ load x[26] from s26
eor r5, r5, r6    @ x[30] = x[30] ^ x[26]
eor r3, r3, r4    @ x[28] = x[28] ^ x[6]
vmov r7, s31    @ load x[31] from s31
eor r3, r3, r7    @ x[28] = x[28] ^ x[31]
vmov r8, s4    @ load x[4] from s4
eor r8, r8, r0    @ x[4] = x[4] ^ x[15]
vmov r9, s22    @ load x[22] from s22
eor r9, r9, r5    @ x[22] = x[22] ^ x[30]
vmov r10, s12    @ load x[12] from s12
vmov r11, s14    @ load x[14] from s14
eor r10, r10, r11    @ x[12] = x[12] ^ x[14]
vmov s12, r10    @ spill x[12] from r10
vmov r10, s29    @ load x[29] from s29
eor r1, r1, r10    @ x[13] = x[13] ^ x[29]
vmov s14, r11    @ spill x[14] from r11
vmov r11, s23    @ load x[23] from s23
vmov s30, r5    @ spill x[30] from r5
vmov r5, s17    @ load x[17] from s17
eor r11, r11, r5    @ x[23] = x[23] ^ x[17]
vmov s5, r2    @ spill x[5] from r2
vmov r2, s16    @ load x[16] from s16
eor r2, r2, r10    @ x[16] = x[16] ^ x[29]
vmov s29, r10    @ spill x[29] from r10
vmov r10, s1    @ load x[1] from s1
eor r10, r10, r4    @ x[1] = x[1] ^ x[6]
vmov s1, r10    @ spill x[1] from r10
vmov r10, s8    @ load x[8] from s8
eor r10, r10, r7    @ x[8] = x[8] ^ x[31]
vmov s31, r7    @ spill x[31] from r7
vmov r7, s21    @ load x[21] from s21
eor r7, r7, r5    @ x[21] = x[21] ^ x[17]
vmov s8, r10    @ spill x[8] from r10
vmov r10, s7    @ load x[7] from s7
eor r10, r10, r3    @ x[7] = x[7] ^ x[28]
vmov s7, r10    @ spill x[7] from r10
vmov r10, s20    @ load x[20] from s20
eor r5, r5, r10    @ x[17] = x[17] ^ x[20]
eor r9, r9, r5    @ x[22] = x[22] ^ x[17]
eor r5, r5, r0    @ x[17] = x[17] ^ x[15]
eor r3, r3, r10    @ x[28] = x[28] ^ x[20]
vmov s17, r5    @ spill x[17] from r5
vmov r5, s18    @ load x[18] from s18
eor r11, r11, r5    @ x[23] = x[23] ^ x[18]
eor r5, r5, r6    @ x[18] = x[18] ^ x[26]
vmov s26, r6    @ spill x[26] from r6
vmov r6, s27    @ load x[27] from s27
eor r2, r2, r6    @ x[16] = x[16] ^ x[27]
eor r0, r0, r7    @ x[15] = x[15] ^ x[21]
vmov s21, r7    @ spill x[21] from r7
vmov r7, s25    @ load x[25] from s25
eor r7, r7, r0    @ x[25] = x[25] ^ x[15]
eor r0, r0, r9    @ x[15] = x[15] ^ x[22]
eor r5, r5, r9    @ x[18] = x[18] ^ x[22]
vmov s25, r7    @ spill x[25] from r7
vmov r7, s10    @ load x[10] from s10
eor r7, r7, r11    @ x[10] = x[10] ^ x[23]
vmov s18, r5    @ spill x[18] from r5
vmov r5, s19    @ load x[19] from s19
eor r5, r5, r9    @ x[19] = x[19] ^ x[22]
eor r6, r6, r9    @ x[27] = x[27] ^ x[22]
eor r1, r1, r10    @ x[13] = x[13] ^ x[20]
vmov s13, r1    @ spill x[13] from r1
vmov r1, s5    @ load x[5] from s5
eor r1, r1, r8    @ x[5] = x[5] ^ x[4]
eor r8, r8, r6    @ x[4] = x[4] ^ x[27]
eor r4, r4, r2    @ x[6] = x[6] ^ x[16]
eor r8, r8, r3    @ x[4] = x[4] ^ x[28]
eor r11, r11, r0    @ x[23] = x[23] ^ x[15]
vmov s23, r11    @ spill x[23] from r11
vmov r11, s24    @ load x[24] from s24
vmov s15, r0    @ spill x[15] from r0
vmov r0, s30    @ load x[30] from s30
eor r11, r11, r0    @ x[24] = x[24] ^ x[30]
vmov s6, r4    @ spill x[6] from r4
vmov r4, s14    @ load x[14] from s14
vmov s28, r3    @ spill x[28] from r3
vmov r3, s9    @ load x[9] from s9
eor r4, r4, r3    @ x[14] = x[14] ^ x[9]
vmov s10, r7    @ spill x[10] from r7
vmov r7, s12    @ load x[12] from s12
eor r3, r3, r7    @ x[9] = x[9] ^ x[12]
vmov s12, r7    @ spill x[12] from r7
vmov r7, s0    @ load x[0] from s0
eor r9, r9, r7    @ x[22] = x[22] ^ x[0]
vmov s22, r9    @ spill x[22] from r9
vmov r9, s21    @ load x[21] from s21
eor r11, r11, r9    @ x[24] = x[24] ^ x[21]
vmov s24, r11    @ spill x[24] from r11
vmov r11, s26    @ load x[26] from s26
eor r6, r6, r11    @ x[27] = x[27] ^ x[26]
eor r11, r11, r10    @ x[26] = x[26] ^ x[20]
eor r10, r10, r3    @ x[20] = x[20] ^ x[9]
eor r7, r7, r9    @ x[0] = x[0] ^ x[21]
eor r4, r4, r1    @ x[14] = x[14] ^ x[5]
eor r0, r0, r2    @ x[30] = x[30] ^ x[16]
vmov s27, r6    @ spill x[27] from r6
vmov r6, s8    @ load x[8] from s8
eor r2, r2, r6    @ x[16] = x[16] ^ x[8]
eor r1, r1, r6    @ x[5] = x[5] ^ x[8]
eor r9, r9, r1    @ x[21] = x[21] ^ x[5]
vmov s30, r0    @ spill x[30] from r0
vmov r0, s7    @ load x[7] from s7
eor r1, r1, r0    @ x[5] = x[5] ^ x[7]
vmov s5, r1    @ spill x[5] from r1
vmov r1, s17    @ load x[17] from s17
eor r7, r7, r1    @ x[0] = x[0] ^ x[17]
eor r3, r3, r8    @ x[9] = x[9] ^ x[4]
eor r8, r8, r5    @ x[4] = x[4] ^ x[19]
vmov s9, r3    @ spill x[9] from r3
vmov r3, s31    @ load x[31] from s31
eor r9, r9, r3    @ x[21] = x[21] ^ x[31]
eor r2, r2, r10    @ x[16] = x[16] ^ x[20]
vmov s16, r2    @ spill x[16] from r2
vmov r2, s1    @ load x[1] from s1
eor r10, r10, r2    @ x[20] = x[20] ^ x[1]
eor r5, r5, r3    @ x[19] = x[19] ^ x[31]
eor r11, r11, r0    @ x[26] = x[26] ^ x[7]
vmov s20, r10    @ spill x[20] from r10
vmov r10, s2    @ load x[2] from s2
eor r2, r2, r10    @ x[1] = x[1] ^ x[2]
eor r2, r2, r0    @ x[1] = x[1] ^ x[7]
vmov s1, r2    @ spill x[1] from r2
vmov r2, s10    @ load x[10] from s10
eor r0, r0, r2    @ x[7] = x[7] ^ x[10]
eor r2, r2, r6    @ x[10] = x[10] ^ x[8]
eor r6, r6, r1    @ x[8] = x[8] ^ x[17]
vmov s8, r6    @ spill x[8] from r6
vmov r6, s22    @ load x[22] from s22
eor r6, r6, r10    @ x[22] = x[22] ^ x[2]
vmov s26, r11    @ spill x[26] from r11
vmov r11, s18    @ load x[18] from s18
eor r10, r10, r11    @ x[2] = x[2] ^ x[18]
vmov s7, r0    @ spill x[7] from r0
vmov r0, s28    @ load x[28] from s28
eor r3, r3, r0    @ x[31] = x[31] ^ x[28]
vmov s19, r5    @ spill x[19] from r5
vmov r5, s13    @ load x[13] from s13
eor r1, r1, r5    @ x[17] = x[17] ^ x[13]
vmov s10, r2    @ spill x[10] from r2
vmov r2, s25    @ load x[25] from s25
eor r3, r3, r2    @ x[31] = x[31] ^ x[25]
vmov s31, r3    @ spill x[31] from r3
vmov r3, s12    @ load x[12] from s12
eor r5, r5, r3    @ x[13] = x[13] ^ x[12]
eor r0, r0, r11    @ x[28] = x[28] ^ x[18]
vmov s13, r5    @ spill x[13] from r5
vmov r5, s6    @ load x[6] from s6
eor r0, r0, r5    @ x[28] = x[28] ^ x[6]
vmov s28, r0    @ spill x[28] from r0
vmov r0, s29    @ load x[29] from s29
eor r7, r7, r0    @ x[0] = x[0] ^ x[29]
vmov s0, r7    @ spill x[0] from r7
vmov r7, s15    @ load x[15] from s15
eor r8, r8, r7    @ x[4] = x[4] ^ x[15]
eor r10, r10, r0    @ x[2] = x[2] ^ x[29]
eor r5, r5, r7    @ x[6] = x[6] ^ x[15]
eor r7, r7, r0    @ x[15] = x[15] ^ x[29]
eor r4, r4, r9    @ x[14] = x[14] ^ x[21]
vmov s21, r9    @ spill x[21] from r9
vmov r9, s11    @ load x[11] from s11
eor r6, r6, r9    @ x[22] = x[22] ^ x[11]
eor r11, r11, r9    @ x[18] = x[18] ^ x[11]
eor r1, r1, r4    @ x[17] = x[17] ^ x[14]
eor r4, r4, r9    @ x[14] = x[14] ^ x[11]
vmov s22, r6    @ spill x[22] from r6
vmov r6, s10    @ load x[10] from s10
eor r1, r1, r6    @ x[17] = x[17] ^ x[10]
vmov s18, r11    @ spill x[18] from r11
vmov r11, s23    @ load x[23] from s23
eor r0, r0, r11    @ x[29] = x[29] ^ x[23]
vmov s17, r1    @ spill x[17] from r1
vmov r1, s30    @ load x[30] from s30
eor r6, r6, r1    @ x[10] = x[10] ^ x[30]
eor r11, r11, r1    @ x[23] = x[23] ^ x[30]
eor r1, r1, r9    @ x[30] = x[30] ^ x[11]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s24    @ load x[24] from s24
eor r3, r3, r6    @ x[12] = x[12] ^ x[24]
eor r9, r9, r2    @ x[11] = x[11] ^ x[25]
vmov s25, r2    @ spill x[25] from r2
vmov r2, s19    @ load x[19] from s19
vmov s24, r6    @ spill x[24] from r6
vmov r6, s3    @ load x[3] from s3
eor r2, r2, r6    @ x[19] = x[19] ^ x[3]
vmov s3, r6    @ spill x[3] from r6
vmov r6, s7    @ load x[7] from s7
eor r10, r10, r6    @ x[2] = x[2] ^ x[7]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s16    @ load x[16] from s16
eor r2, r2, r4    @ x[16] = x[16] ^ x[14]
vmov s16, r2    @ spill x[16] from r2
vmov r2, s31    @ load x[31] from s31
eor r4, r4, r2    @ x[14] = x[14] ^ x[31]
eor r2, r2, r7    @ x[31] = x[31] ^ x[15]
eor r9, r9, r0    @ x[11] = x[11] ^ x[29]
eor r1, r1, r5    @ x[30] = x[30] ^ x[6]
vmov s14, r4    @ spill x[14] from r4
vmov r4, s26    @ load x[26] from s26
eor r7, r7, r4    @ x[15] = x[15] ^ x[26]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s1    @ load x[1] from s1
eor r7, r7, r2    @ x[15] = x[15] ^ x[1]
eor r8, r8, r0    @ x[4] = x[4] ^ x[29]
eor r6, r6, r11    @ x[7] = x[7] ^ x[23]
eor r5, r5, r6    @ x[6] = x[6] ^ x[7]
vmov s4, r8    @ spill x[4] from r8
vmov r8, s9    @ load x[9] from s9
eor r6, r6, r8    @ x[7] = x[7] ^ x[9]
vmov s9, r8    @ spill x[9] from r8
vmov r8, s28    @ load x[28] from s28
eor r8, r8, r0    @ x[28] = x[28] ^ x[29]
eor r5, r5, r3    @ x[6] = x[6] ^ x[12]
vmov s6, r5    @ spill x[6] from r5
vmov r5, s13    @ load x[13] from s13
eor r10, r10, r5    @ x[2] = x[2] ^ x[13]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s5    @ load x[5] from s5
vmov s13, r5    @ spill x[13] from r5
vmov r5, s0    @ load x[0] from s0
eor r10, r10, r5    @ x[5] = x[5] ^ x[0]
eor r1, r1, r7    @ x[30] = x[30] ^ x[15]
vmov s5, r10    @ spill x[5] from r10
vmov r10, s19    @ load x[19] from s19
eor r11, r11, r10    @ x[23] = x[23] ^ x[19]
vmov s30, r1    @ spill x[30] from r1
vmov r1, s16    @ load x[16] from s16
eor r6, r6, r1    @ x[7] = x[7] ^ x[16]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s17    @ load x[17] from s17
eor r6, r6, r8    @ x[17] = x[17] ^ x[28]
vmov s17, r6    @ spill x[17] from r6
vmov r6, s27    @ load x[27] from s27
eor r9, r9, r6    @ x[11] = x[11] ^ x[27]
eor r7, r7, r10    @ x[15] = x[15] ^ x[19]
vmov s11, r9    @ spill x[11] from r9
vmov r9, s3    @ load x[3] from s3
eor r4, r4, r9    @ x[26] = x[26] ^ x[3]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s24    @ load x[24] from s24
eor r0, r0, r7    @ x[29] = x[29] ^ x[24]
eor r10, r10, r7    @ x[19] = x[19] ^ x[24]
vmov s29, r0    @ spill x[29] from r0
vmov r0, s10    @ load x[10] from s10
vmov s26, r4    @ spill x[26] from r4
vmov r4, s18    @ load x[18] from s18
eor r0, r0, r4    @ x[10] = x[10] ^ x[18]
eor r10, r10, r4    @ x[19] = x[19] ^ x[18]
eor r9, r9, r6    @ x[3] = x[3] ^ x[27]
eor r6, r6, r2    @ x[27] = x[27] ^ x[1]
vmov s10, r0    @ spill x[10] from r0
vmov r0, s21    @ load x[21] from s21
eor r7, r7, r0    @ x[24] = x[24] ^ x[21]
eor r2, r2, r0    @ x[1] = x[1] ^ x[21]
vmov s27, r6    @ spill x[27] from r6
vmov r6, s22    @ load x[22] from s22
eor r0, r0, r6    @ x[21] = x[21] ^ x[22]
eor r0, r0, r3    @ x[21] = x[21] ^ x[12]
vmov s21, r0    @ spill x[21] from r0
vmov r0, s25    @ load x[25] from s25
eor r2, r2, r0    @ x[1] = x[1] ^ x[25]
eor r3, r3, r8    @ x[12] = x[12] ^ x[28]
eor r8, r8, r5    @ x[28] = x[28] ^ x[0]
eor r0, r0, r6    @ x[25] = x[25] ^ x[22]
eor r6, r6, r4    @ x[22] = x[22] ^ x[18]
eor r4, r4, r1    @ x[18] = x[18] ^ x[16]
eor r6, r6, r11    @ x[22] = x[22] ^ x[23]
vmov s25, r0    @ spill x[25] from r0
vmov r0, s13    @ load x[13] from s13
eor r4, r4, r0    @ x[18] = x[18] ^ x[13]
vmov s22, r6    @ spill x[22] from r6
vmov r6, s9    @ load x[9] from s9
eor r5, r5, r6    @ x[0] = x[0] ^ x[9]
eor r11, r11, r0    @ x[23] = x[23] ^ x[13]
eor r3, r3, r1    @ x[12] = x[12] ^ x[16]
vmov s12, r3    @ spill x[12] from r3
vmov r3, s31    @ load x[31] from s31
eor r0, r0, r3    @ x[13] = x[13] ^ x[31]
eor r3, r3, r6    @ x[31] = x[31] ^ x[9]
vmov s31, r3    @ spill x[31] from r3
vmov r3, s6    @ load x[6] from s6
eor r1, r1, r3    @ x[16] = x[16] ^ x[6]
eor r2, r2, r1    @ x[1] = x[1] ^ x[16]
vmov s1, r2    @ spill x[1] from r2
vmov r2, s14    @ load x[14] from s14
eor r6, r6, r2    @ x[9] = x[9] ^ x[14]
vmov s6, r3    @ spill x[6] from r3
vmov r3, s30    @ load x[30] from s30
eor r10, r10, r3    @ x[19] = x[19] ^ x[30]
vmov s19, r10    @ spill x[19] from r10
vmov r10, s26    @ load x[26] from s26
eor r5, r5, r10    @ x[0] = x[0] ^ x[26]
eor r10, r10, r3    @ x[26] = x[26] ^ x[30]
eor r1, r1, r2    @ x[16] = x[16] ^ x[14]
eor r2, r2, r3    @ x[14] = x[14] ^ x[30]
eor r11, r11, r9    @ x[23] = x[23] ^ x[3]
vmov s14, r2    @ spill x[14] from r2
vmov r2, s8    @ load x[8] from s8
eor r9, r9, r2    @ x[3] = x[3] ^ x[8]
eor r6, r6, r0    @ x[9] = x[9] ^ x[13]
eor r2, r2, r7    @ x[8] = x[8] ^ x[24]
vmov s24, r7    @ spill x[24] from r7
vmov r7, s15    @ load x[15] from s15
eor r7, r7, r8    @ x[15] = x[15] ^ x[28]
vmov s8, r2    @ spill x[8] from r2
vmov r2, s11    @ load x[11] from s11
eor r3, r3, r2    @ x[30] = x[30] ^ x[11]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s27    @ load x[27] from s27
eor r7, r7, r3    @ x[27] = x[27] ^ x[30]
vmov s16, r1    @ spill x[16] from r1
vmov r1, s21    @ load x[21] from s21
eor r1, r1, r0    @ x[21] = x[21] ^ x[13]
eor r0, r0, r11    @ x[13] = x[13] ^ x[23]
vmov s13, r0    @ spill x[13] from r0
vmov r0, s17    @ load x[17] from s17
eor r0, r0, r4    @ x[17] = x[17] ^ x[18]
vmov s18, r4    @ spill x[18] from r4
vmov r4, s22    @ load x[22] from s22
vmov s27, r7    @ spill x[27] from r7
vmov r7, s29    @ load x[29] from s29
eor r4, r4, r7    @ x[22] = x[22] ^ x[29]
eor r5, r5, r0    @ x[0] = x[0] ^ x[17]
eor r9, r9, r6    @ x[3] = x[3] ^ x[9]
vmov s0, r5    @ spill x[0] from r5
vmov r5, s12    @ load x[12] from s12
eor r5, r5, r2    @ x[12] = x[12] ^ x[11]
eor r8, r8, r3    @ x[28] = x[28] ^ x[30]
eor r7, r7, r11    @ x[29] = x[29] ^ x[23]
eor r1, r1, r10    @ x[21] = x[21] ^ x[26]
vmov s23, r11    @ spill x[23] from r11
vmov r11, s27    @ load x[27] from s27
vmov s28, r8    @ spill x[28] from r8
vmov r8, s13    @ load x[13] from s13
eor r11, r11, r8    @ x[27] = x[27] ^ x[13]
vmov s3, r9    @ spill x[3] from r9
vmov r9, s16    @ load x[16] from s16
eor r9, r9, r11    @ x[16] = x[16] ^ x[27]
eor r2, r2, r0    @ x[11] = x[11] ^ x[17]
vmov s9, r6    @ spill x[9] from r6
vmov r6, s14    @ load x[14] from s14
eor r6, r6, r1    @ x[14] = x[14] ^ x[21]
vmov s30, r3    @ spill x[30] from r3
vmov r3, s15    @ load x[15] from s15
vmov s26, r10    @ spill x[26] from r10
vmov r10, s6    @ load x[6] from s6
eor r3, r3, r10    @ x[15] = x[15] ^ x[6]
vmov s11, r2    @ spill x[11] from r2
vmov r2, s1    @ load x[1] from s1
vmov s17, r0    @ spill x[17] from r0
vmov r0, s31    @ load x[31] from s31
eor r2, r2, r0    @ x[1] = x[1] ^ x[31]
vmov s22, r4    @ spill x[22] from r4
vmov r4, s8    @ load x[8] from s8
eor r4, r4, r11    @ x[8] = x[8] ^ x[27]
eor r10, r10, r1    @ x[6] = x[6] ^ x[21]
vmov s21, r1    @ spill x[21] from r1
vmov r1, s25    @ load x[25] from s25
vmov s29, r7    @ spill x[29] from r7
vmov r7, s4    @ load x[4] from s4
eor r1, r1, r7    @ x[25] = x[25] ^ x[4]
vstr.32  s9, [r14, #0]    @ y[0] = x[9] from s9
vstr.32  s7, [r14, #4]    @ y[1] = x[7] from s7
vstr.32  s24, [r14, #8]    @ y[2] = x[24] from s24
str  r0, [r14, #12]    @ y[3] = x[31] from r0
vstr.32  s21, [r14, #16]    @ y[4] = x[21] from s21
str  r6, [r14, #20]    @ y[5] = x[14] from r6
vstr.32  s0, [r14, #24]    @ y[6] = x[0] from s0
vstr.32  s26, [r14, #28]    @ y[7] = x[26] from s26
vstr.32  s5, [r14, #32]    @ y[8] = x[5] from s5
vstr.32  s20, [r14, #36]    @ y[9] = x[20] from s20
vstr.32  s2, [r14, #40]    @ y[10] = x[2] from s2
vstr.32  s3, [r14, #44]    @ y[11] = x[3] from s3
str  r3, [r14, #48]    @ y[12] = x[15] from r3
str  r1, [r14, #52]    @ y[13] = x[25] from r1
str  r2, [r14, #56]    @ y[14] = x[1] from r2
str  r10, [r14, #60]    @ y[15] = x[6] from r10
vstr.32  s22, [r14, #64]    @ y[16] = x[22] from s22
vstr.32  s30, [r14, #68]    @ y[17] = x[30] from s30
vstr.32  s19, [r14, #72]    @ y[18] = x[19] from s19
vstr.32  s28, [r14, #76]    @ y[19] = x[28] from s28
str  r4, [r14, #80]    @ y[20] = x[8] from r4
str  r8, [r14, #84]    @ y[21] = x[13] from r8
str  r11, [r14, #88]    @ y[22] = x[27] from r11
str  r7, [r14, #92]    @ y[23] = x[4] from r7
vstr.32  s10, [r14, #96]    @ y[24] = x[10] from s10
vstr.32  s17, [r14, #100]    @ y[25] = x[17] from s17
vstr.32  s11, [r14, #104]    @ y[26] = x[11] from s11
str  r9, [r14, #108]    @ y[27] = x[16] from r9
vstr.32  s23, [r14, #112]    @ y[28] = x[23] from s23
vstr.32  s18, [r14, #116]    @ y[29] = x[18] from s18
str  r5, [r14, #120]    @ y[30] = x[12] from r5
vstr.32  s29, [r14, #124]    @ y[31] = x[29] from s29
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v19, .-gft_mul_v19
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v20
.type gft_mul_v20, %function
.align 2
gft_mul_v20:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s10    @ load x[10] from s10
vmov r1, s26    @ load x[26] from s26
eor r0, r0, r1    @ x[10] = x[10] ^ x[26]
vmov r2, s31    @ load x[31] from s31
eor r1, r1, r2    @ x[26] = x[26] ^ x[31]
vmov r3, s20    @ load x[20] from s20
vmov r4, s18    @ load x[18] from s18
eor r3, r3, r4    @ x[20] = x[20] ^ x[18]
vmov r5, s27    @ load x[27] from s27
vmov r6, s24    @ load x[24] from s24
eor r5, r5, r6    @ x[27] = x[27] ^ x[24]
vmov r7, s5    @ load x[5] from s5
eor r4, r4, r7    @ x[18] = x[18] ^ x[5]
vmov r8, s21    @ load x[21] from s21
eor r3, r3, r8    @ x[20] = x[20] ^ x[21]
vmov r9, s30    @ load x[30] from s30
eor r4, r4, r9    @ x[18] = x[18] ^ x[30]
eor r0, r0, r4    @ x[10] = x[10] ^ x[18]
vmov r10, s16    @ load x[16] from s16
eor r4, r4, r10    @ x[18] = x[18] ^ x[16]
vmov r11, s23    @ load x[23] from s23
vmov s21, r8    @ spill x[21] from r8
vmov r8, s28    @ load x[28] from s28
eor r11, r11, r8    @ x[23] = x[23] ^ x[28]
vmov s28, r8    @ spill x[28] from r8
vmov r8, s29    @ load x[29] from s29
eor r8, r8, r5    @ x[29] = x[29] ^ x[27]
vmov s27, r5    @ spill x[27] from r5
vmov r5, s15    @ load x[15] from s15
eor r5, r5, r6    @ x[15] = x[15] ^ x[24]
vmov s16, r10    @ spill x[16] from r10
vmov r10, s22    @ load x[22] from s22
eor r10, r10, r4    @ x[22] = x[22] ^ x[18]
vmov s22, r10    @ spill x[22] from r10
vmov r10, s2    @ load x[2] from s2
vmov s26, r1    @ spill x[26] from r1
vmov r1, s6    @ load x[6] from s6
eor r10, r10, r1    @ x[2] = x[2] ^ x[6]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s13    @ load x[13] from s13
vmov s10, r0    @ spill x[10] from r0
vmov r0, s25    @ load x[25] from s25
eor r10, r10, r0    @ x[13] = x[13] ^ x[25]
eor r9, r9, r6    @ x[30] = x[30] ^ x[24]
vmov s25, r0    @ spill x[25] from r0
vmov r0, s7    @ load x[7] from s7
eor r7, r7, r0    @ x[5] = x[5] ^ x[7]
vmov s7, r0    @ spill x[7] from r0
vmov r0, s4    @ load x[4] from s4
eor r0, r0, r1    @ x[4] = x[4] ^ x[6]
eor r1, r1, r3    @ x[6] = x[6] ^ x[20]
eor r6, r6, r0    @ x[24] = x[24] ^ x[4]
vmov s13, r10    @ spill x[13] from r10
vmov r10, s1    @ load x[1] from s1
vmov s5, r7    @ spill x[5] from r7
vmov r7, s14    @ load x[14] from s14
eor r10, r10, r7    @ x[1] = x[1] ^ x[14]
eor r7, r7, r11    @ x[14] = x[14] ^ x[23]
eor r0, r0, r2    @ x[4] = x[4] ^ x[31]
eor r6, r6, r4    @ x[24] = x[24] ^ x[18]
vmov s4, r0    @ spill x[4] from r0
vmov r0, s17    @ load x[17] from s17
vmov s23, r11    @ spill x[23] from r11
vmov r11, s0    @ load x[0] from s0
eor r0, r0, r11    @ x[17] = x[17] ^ x[0]
vmov s18, r4    @ spill x[18] from r4
vmov r4, s3    @ load x[3] from s3
eor r4, r4, r1    @ x[3] = x[3] ^ x[6]
eor r0, r0, r3    @ x[17] = x[17] ^ x[20]
eor r3, r3, r9    @ x[20] = x[20] ^ x[30]
vmov s3, r4    @ spill x[3] from r4
vmov r4, s9    @ load x[9] from s9
eor r4, r4, r5    @ x[9] = x[9] ^ x[15]
vmov s9, r4    @ spill x[9] from r4
vmov r4, s10    @ load x[10] from s10
eor r4, r4, r8    @ x[10] = x[10] ^ x[29]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s5    @ load x[5] from s5
vmov s17, r0    @ spill x[17] from r0
vmov r0, s26    @ load x[26] from s26
eor r1, r1, r0    @ x[5] = x[5] ^ x[26]
eor r4, r4, r5    @ x[10] = x[10] ^ x[15]
vmov s10, r4    @ spill x[10] from r4
vmov r4, s19    @ load x[19] from s19
eor r8, r8, r4    @ x[29] = x[29] ^ x[19]
vmov s15, r5    @ spill x[15] from r5
vmov r5, s16    @ load x[16] from s16
eor r4, r4, r5    @ x[19] = x[19] ^ x[16]
vmov s20, r3    @ spill x[20] from r3
vmov r3, s22    @ load x[22] from s22
eor r8, r8, r3    @ x[29] = x[29] ^ x[22]
vmov s29, r8    @ spill x[29] from r8
vmov r8, s21    @ load x[21] from s21
eor r4, r4, r8    @ x[19] = x[19] ^ x[21]
eor r5, r5, r6    @ x[16] = x[16] ^ x[24]
vmov s24, r6    @ spill x[24] from r6
vmov r6, s28    @ load x[28] from s28
eor r4, r4, r6    @ x[19] = x[19] ^ x[28]
eor r5, r5, r9    @ x[16] = x[16] ^ x[30]
eor r6, r6, r10    @ x[28] = x[28] ^ x[1]
eor r10, r10, r3    @ x[1] = x[1] ^ x[22]
vmov s28, r6    @ spill x[28] from r6
vmov r6, s13    @ load x[13] from s13
eor r9, r9, r6    @ x[30] = x[30] ^ x[13]
vmov s19, r4    @ spill x[19] from r4
vmov r4, s8    @ load x[8] from s8
eor r9, r9, r4    @ x[30] = x[30] ^ x[8]
eor r3, r3, r0    @ x[22] = x[22] ^ x[26]
vmov s16, r5    @ spill x[16] from r5
vmov r5, s2    @ load x[2] from s2
eor r10, r10, r5    @ x[1] = x[1] ^ x[2]
eor r10, r10, r2    @ x[1] = x[1] ^ x[31]
eor r2, r2, r8    @ x[31] = x[31] ^ x[21]
eor r8, r8, r4    @ x[21] = x[21] ^ x[8]
eor r4, r4, r7    @ x[8] = x[8] ^ x[14]
vmov s21, r8    @ spill x[21] from r8
vmov r8, s11    @ load x[11] from s11
eor r7, r7, r8    @ x[14] = x[14] ^ x[11]
vmov s8, r4    @ spill x[8] from r4
vmov r4, s18    @ load x[18] from s18
eor r4, r4, r11    @ x[18] = x[18] ^ x[0]
eor r8, r8, r11    @ x[11] = x[11] ^ x[0]
eor r11, r11, r6    @ x[0] = x[0] ^ x[13]
eor r6, r6, r4    @ x[13] = x[13] ^ x[18]
eor r6, r6, r8    @ x[13] = x[13] ^ x[11]
vmov s18, r4    @ spill x[18] from r4
vmov r4, s12    @ load x[12] from s12
eor r11, r11, r4    @ x[0] = x[0] ^ x[12]
eor r4, r4, r5    @ x[12] = x[12] ^ x[2]
vmov s0, r11    @ spill x[0] from r11
vmov r11, s25    @ load x[25] from s25
eor r7, r7, r11    @ x[14] = x[14] ^ x[25]
eor r11, r11, r5    @ x[25] = x[25] ^ x[2]
eor r5, r5, r1    @ x[2] = x[2] ^ x[5]
eor r8, r8, r1    @ x[11] = x[11] ^ x[5]
vmov s14, r7    @ spill x[14] from r7
vmov r7, s27    @ load x[27] from s27
eor r1, r1, r7    @ x[5] = x[5] ^ x[27]
vmov s5, r1    @ spill x[5] from r1
vmov r1, s24    @ load x[24] from s24
eor r2, r2, r1    @ x[31] = x[31] ^ x[24]
vmov s25, r11    @ spill x[25] from r11
vmov r11, s20    @ load x[20] from s20
eor r0, r0, r11    @ x[26] = x[26] ^ x[20]
eor r5, r5, r0    @ x[2] = x[2] ^ x[26]
eor r4, r4, r7    @ x[12] = x[12] ^ x[27]
eor r9, r9, r7    @ x[30] = x[30] ^ x[27]
vmov s2, r5    @ spill x[2] from r5
vmov r5, s17    @ load x[17] from s17
eor r5, r5, r3    @ x[17] = x[17] ^ x[22]
vmov s12, r4    @ spill x[12] from r4
vmov r4, s23    @ load x[23] from s23
eor r0, r0, r4    @ x[26] = x[26] ^ x[23]
vmov s26, r0    @ spill x[26] from r0
vmov r0, s6    @ load x[6] from s6
eor r7, r7, r0    @ x[27] = x[27] ^ x[6]
eor r7, r7, r1    @ x[27] = x[27] ^ x[24]
eor r2, r2, r6    @ x[31] = x[31] ^ x[13]
vmov s13, r6    @ spill x[13] from r6
vmov r6, s15    @ load x[15] from s15
eor r11, r11, r6    @ x[20] = x[20] ^ x[15]
vmov s15, r6    @ spill x[15] from r6
vmov r6, s29    @ load x[29] from s29
eor r0, r0, r6    @ x[6] = x[6] ^ x[29]
eor r1, r1, r4    @ x[24] = x[24] ^ x[23]
eor r4, r4, r8    @ x[23] = x[23] ^ x[11]
eor r4, r4, r5    @ x[23] = x[23] ^ x[17]
vmov s17, r5    @ spill x[17] from r5
vmov r5, s4    @ load x[4] from s4
eor r1, r1, r5    @ x[24] = x[24] ^ x[4]
vmov s29, r6    @ spill x[29] from r6
vmov r6, s10    @ load x[10] from s10
eor r1, r1, r6    @ x[24] = x[24] ^ x[10]
eor r4, r4, r10    @ x[23] = x[23] ^ x[1]
eor r4, r4, r9    @ x[23] = x[23] ^ x[30]
vmov s1, r10    @ spill x[1] from r10
vmov r10, s9    @ load x[9] from s9
eor r4, r4, r10    @ x[23] = x[23] ^ x[9]
eor r4, r4, r2    @ x[23] = x[23] ^ x[31]
eor r5, r5, r3    @ x[4] = x[4] ^ x[22]
vmov s23, r4    @ spill x[23] from r4
vmov r4, s7    @ load x[7] from s7
eor r3, r3, r4    @ x[22] = x[22] ^ x[7]
vmov s24, r1    @ spill x[24] from r1
vmov r1, s25    @ load x[25] from s25
eor r3, r3, r1    @ x[22] = x[22] ^ x[25]
vmov s11, r8    @ spill x[11] from r8
vmov r8, s16    @ load x[16] from s16
eor r1, r1, r8    @ x[25] = x[25] ^ x[16]
vmov s7, r4    @ spill x[7] from r4
vmov r4, s3    @ load x[3] from s3
eor r2, r2, r4    @ x[31] = x[31] ^ x[3]
eor r1, r1, r7    @ x[25] = x[25] ^ x[27]
eor r3, r3, r0    @ x[22] = x[22] ^ x[6]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s26    @ load x[26] from s26
eor r0, r0, r2    @ x[6] = x[6] ^ x[26]
vmov s22, r3    @ spill x[22] from r3
vmov r3, s8    @ load x[8] from s8
eor r7, r7, r3    @ x[27] = x[27] ^ x[8]
vmov s27, r7    @ spill x[27] from r7
vmov r7, s18    @ load x[18] from s18
eor r0, r0, r7    @ x[6] = x[6] ^ x[18]
eor r1, r1, r9    @ x[25] = x[25] ^ x[30]
eor r1, r1, r4    @ x[25] = x[25] ^ x[3]
eor r5, r5, r10    @ x[4] = x[4] ^ x[9]
vmov s4, r5    @ spill x[4] from r5
vmov r5, s19    @ load x[19] from s19
eor r2, r2, r5    @ x[26] = x[26] ^ x[19]
eor r4, r4, r5    @ x[3] = x[3] ^ x[19]
eor r5, r5, r9    @ x[19] = x[19] ^ x[30]
vmov s19, r5    @ spill x[19] from r5
vmov r5, s0    @ load x[0] from s0
eor r9, r9, r5    @ x[30] = x[30] ^ x[0]
eor r5, r5, r10    @ x[0] = x[0] ^ x[9]
eor r10, r10, r2    @ x[9] = x[9] ^ x[26]
eor r10, r10, r4    @ x[9] = x[9] ^ x[3]
vmov s6, r0    @ spill x[6] from r0
vmov r0, s5    @ load x[5] from s5
eor r9, r9, r0    @ x[30] = x[30] ^ x[5]
eor r0, r0, r3    @ x[5] = x[5] ^ x[8]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s21    @ load x[21] from s21
eor r10, r10, r9    @ x[9] = x[9] ^ x[21]
eor r10, r10, r11    @ x[9] = x[9] ^ x[20]
eor r4, r4, r7    @ x[3] = x[3] ^ x[18]
eor r3, r3, r2    @ x[8] = x[8] ^ x[26]
eor r7, r7, r11    @ x[18] = x[18] ^ x[20]
eor r2, r2, r8    @ x[26] = x[26] ^ x[16]
eor r8, r8, r6    @ x[16] = x[16] ^ x[10]
vmov s9, r10    @ spill x[9] from r10
vmov r10, s13    @ load x[13] from s13
eor r6, r6, r10    @ x[10] = x[10] ^ x[13]
vmov s25, r1    @ spill x[25] from r1
vmov r1, s7    @ load x[7] from s7
eor r11, r11, r1    @ x[20] = x[20] ^ x[7]
vmov s20, r11    @ spill x[20] from r11
vmov r11, s27    @ load x[27] from s27
vmov s8, r3    @ spill x[8] from r3
vmov r3, s28    @ load x[28] from s28
eor r11, r11, r3    @ x[27] = x[27] ^ x[28]
vmov s26, r2    @ spill x[26] from r2
vmov r2, s11    @ load x[11] from s11
eor r1, r1, r2    @ x[7] = x[7] ^ x[11]
eor r1, r1, r0    @ x[7] = x[7] ^ x[5]
vmov s18, r7    @ spill x[18] from r7
vmov r7, s24    @ load x[24] from s24
eor r1, r1, r7    @ x[7] = x[7] ^ x[24]
eor r2, r2, r3    @ x[11] = x[11] ^ x[28]
vmov s11, r2    @ spill x[11] from r2
vmov r2, s29    @ load x[29] from s29
eor r10, r10, r2    @ x[13] = x[13] ^ x[29]
vmov s13, r10    @ spill x[13] from r10
vmov r10, s15    @ load x[15] from s15
eor r2, r2, r10    @ x[29] = x[29] ^ x[15]
eor r0, r0, r11    @ x[5] = x[5] ^ x[27]
eor r3, r3, r8    @ x[28] = x[28] ^ x[16]
vmov s28, r3    @ spill x[28] from r3
vmov r3, s12    @ load x[12] from s12
eor r8, r8, r3    @ x[16] = x[16] ^ x[12]
eor r0, r0, r4    @ x[5] = x[5] ^ x[3]
vmov s16, r8    @ spill x[16] from r8
vmov r8, s2    @ load x[2] from s2
eor r4, r4, r8    @ x[3] = x[3] ^ x[2]
eor r11, r11, r9    @ x[27] = x[27] ^ x[21]
eor r3, r3, r1    @ x[12] = x[12] ^ x[7]
eor r6, r6, r5    @ x[10] = x[10] ^ x[0]
vmov s7, r1    @ spill x[7] from r1
vmov r1, s17    @ load x[17] from s17
eor r5, r5, r1    @ x[0] = x[0] ^ x[17]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s14    @ load x[14] from s14
eor r6, r6, r1    @ x[14] = x[14] ^ x[17]
vmov s3, r4    @ spill x[3] from r4
vmov r4, s13    @ load x[13] from s13
eor r4, r4, r7    @ x[13] = x[13] ^ x[24]
vmov s24, r7    @ spill x[24] from r7
vmov r7, s18    @ load x[18] from s18
eor r6, r6, r7    @ x[14] = x[14] ^ x[18]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s26    @ load x[26] from s26
eor r7, r7, r6    @ x[18] = x[18] ^ x[26]
vmov s13, r4    @ spill x[13] from r4
vmov r4, s8    @ load x[8] from s8
eor r9, r9, r4    @ x[21] = x[21] ^ x[8]
vmov s17, r1    @ spill x[17] from r1
vmov r1, s25    @ load x[25] from s25
eor r7, r7, r1    @ x[18] = x[18] ^ x[25]
eor r8, r8, r9    @ x[2] = x[2] ^ x[21]
eor r0, r0, r3    @ x[5] = x[5] ^ x[12]
eor r3, r3, r1    @ x[12] = x[12] ^ x[25]
eor r2, r2, r5    @ x[29] = x[29] ^ x[0]
eor r10, r10, r4    @ x[15] = x[15] ^ x[8]
eor r7, r7, r10    @ x[18] = x[18] ^ x[15]
vmov s29, r2    @ spill x[29] from r2
vmov r2, s22    @ load x[22] from s22
eor r1, r1, r2    @ x[25] = x[25] ^ x[22]
eor r10, r10, r8    @ x[15] = x[15] ^ x[2]
vmov s12, r3    @ spill x[12] from r3
vmov r3, s16    @ load x[16] from s16
eor r3, r3, r9    @ x[16] = x[16] ^ x[21]
vmov s5, r0    @ spill x[5] from r0
vmov r0, s20    @ load x[20] from s20
eor r5, r5, r0    @ x[0] = x[0] ^ x[20]
vmov s18, r7    @ spill x[18] from r7
vmov r7, s1    @ load x[1] from s1
eor r8, r8, r7    @ x[2] = x[2] ^ x[1]
eor r6, r6, r2    @ x[26] = x[26] ^ x[22]
vmov s26, r6    @ spill x[26] from r6
vmov r6, s6    @ load x[6] from s6
eor r8, r8, r6    @ x[2] = x[2] ^ x[6]
vmov s6, r6    @ spill x[6] from r6
vmov r6, s4    @ load x[4] from s4
eor r7, r7, r6    @ x[1] = x[1] ^ x[4]
eor r0, r0, r11    @ x[20] = x[20] ^ x[27]
vmov s20, r0    @ spill x[20] from r0
vmov r0, s11    @ load x[11] from s11
eor r9, r9, r0    @ x[21] = x[21] ^ x[11]
eor r5, r5, r4    @ x[0] = x[0] ^ x[8]
vmov s8, r4    @ spill x[8] from r4
vmov r4, s23    @ load x[23] from s23
eor r4, r4, r2    @ x[23] = x[23] ^ x[22]
eor r4, r4, r11    @ x[23] = x[23] ^ x[27]
vmov s27, r11    @ spill x[27] from r11
vmov r11, s3    @ load x[3] from s3
vmov s22, r2    @ spill x[22] from r2
vmov r2, s17    @ load x[17] from s17
eor r11, r11, r2    @ x[3] = x[3] ^ x[17]
eor r9, r9, r10    @ x[21] = x[21] ^ x[15]
vmov s11, r0    @ spill x[11] from r0
vmov r0, s19    @ load x[19] from s19
eor r9, r9, r0    @ x[21] = x[21] ^ x[19]
eor r3, r3, r1    @ x[16] = x[16] ^ x[25]
eor r1, r1, r0    @ x[25] = x[25] ^ x[19]
eor r10, r10, r7    @ x[15] = x[15] ^ x[1]
vmov s25, r1    @ spill x[25] from r1
vmov r1, s28    @ load x[28] from s28
eor r10, r10, r1    @ x[15] = x[15] ^ x[28]
vmov s16, r3    @ spill x[16] from r3
vmov r3, s30    @ load x[30] from s30
eor r3, r3, r1    @ x[30] = x[30] ^ x[28]
vmov s23, r4    @ spill x[23] from r4
vmov r4, s31    @ load x[31] from s31
eor r4, r4, r10    @ x[31] = x[31] ^ x[15]
eor r5, r5, r10    @ x[0] = x[0] ^ x[15]
vmov s3, r11    @ spill x[3] from r11
vmov r11, s18    @ load x[18] from s18
eor r11, r11, r5    @ x[18] = x[18] ^ x[0]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s5    @ load x[5] from s5
eor r6, r6, r2    @ x[4] = x[4] ^ x[5]
eor r7, r7, r9    @ x[1] = x[1] ^ x[21]
vmov s21, r9    @ spill x[21] from r9
vmov r9, s10    @ load x[10] from s10
eor r9, r9, r8    @ x[10] = x[10] ^ x[2]
vmov s2, r8    @ spill x[2] from r8
vmov r8, s12    @ load x[12] from s12
eor r10, r10, r8    @ x[15] = x[15] ^ x[12]
vmov s15, r10    @ spill x[15] from r10
vmov r10, s7    @ load x[7] from s7
eor r10, r10, r5    @ x[7] = x[7] ^ x[0]
vmov s0, r5    @ spill x[0] from r5
vmov r5, s20    @ load x[20] from s20
eor r5, r5, r8    @ x[20] = x[20] ^ x[12]
vmov s1, r7    @ spill x[1] from r7
vmov r7, s9    @ load x[9] from s9
eor r7, r7, r10    @ x[9] = x[9] ^ x[7]
vstr.32  s21, [r14, #0]    @ y[0] = x[21] from s21
str  r8, [r14, #4]    @ y[1] = x[12] from r8
str  r3, [r14, #8]    @ y[2] = x[30] from r3
vstr.32  s1, [r14, #12]    @ y[3] = x[1] from s1
str  r1, [r14, #16]    @ y[4] = x[28] from r1
vstr.32  s11, [r14, #20]    @ y[5] = x[11] from s11
str  r4, [r14, #24]    @ y[6] = x[31] from r4
vstr.32  s16, [r14, #28]    @ y[7] = x[16] from s16
str  r5, [r14, #32]    @ y[8] = x[20] from r5
str  r11, [r14, #36]    @ y[9] = x[18] from r11
str  r7, [r14, #40]    @ y[10] = x[9] from r7
str  r6, [r14, #44]    @ y[11] = x[4] from r6
vstr.32  s15, [r14, #48]    @ y[12] = x[15] from s15
str  r2, [r14, #52]    @ y[13] = x[5] from r2
str  r10, [r14, #56]    @ y[14] = x[7] from r10
str  r9, [r14, #60]    @ y[15] = x[10] from r9
vstr.32  s2, [r14, #64]    @ y[16] = x[2] from s2
vstr.32  s6, [r14, #68]    @ y[17] = x[6] from s6
vstr.32  s3, [r14, #72]    @ y[18] = x[3] from s3
vstr.32  s26, [r14, #76]    @ y[19] = x[26] from s26
vstr.32  s17, [r14, #80]    @ y[20] = x[17] from s17
vstr.32  s27, [r14, #84]    @ y[21] = x[27] from s27
vstr.32  s22, [r14, #88]    @ y[22] = x[22] from s22
vstr.32  s25, [r14, #92]    @ y[23] = x[25] from s25
vstr.32  s23, [r14, #96]    @ y[24] = x[23] from s23
vstr.32  s8, [r14, #100]    @ y[25] = x[8] from s8
vstr.32  s13, [r14, #104]    @ y[26] = x[13] from s13
vstr.32  s0, [r14, #108]    @ y[27] = x[0] from s0
str  r0, [r14, #112]    @ y[28] = x[19] from r0
vstr.32  s29, [r14, #116]    @ y[29] = x[29] from s29
vstr.32  s24, [r14, #120]    @ y[30] = x[24] from s24
vstr.32  s14, [r14, #124]    @ y[31] = x[14] from s14
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v20, .-gft_mul_v20
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v21
.type gft_mul_v21, %function
.align 2
gft_mul_v21:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s19    @ load x[19] from s19
vmov r1, s27    @ load x[27] from s27
eor r0, r0, r1    @ x[19] = x[19] ^ x[27]
vmov r2, s2    @ load x[2] from s2
vmov r3, s31    @ load x[31] from s31
eor r2, r2, r3    @ x[2] = x[2] ^ x[31]
vmov r4, s1    @ load x[1] from s1
vmov r5, s7    @ load x[7] from s7
eor r4, r4, r5    @ x[1] = x[1] ^ x[7]
vmov r6, s26    @ load x[26] from s26
vmov r7, s3    @ load x[3] from s3
eor r6, r6, r7    @ x[26] = x[26] ^ x[3]
vmov r8, s30    @ load x[30] from s30
vmov r9, s18    @ load x[18] from s18
eor r8, r8, r9    @ x[30] = x[30] ^ x[18]
vmov r10, s13    @ load x[13] from s13
eor r10, r10, r4    @ x[13] = x[13] ^ x[1]
vmov r11, s14    @ load x[14] from s14
eor r11, r11, r10    @ x[14] = x[14] ^ x[13]
vmov s27, r1    @ spill x[27] from r1
vmov r1, s6    @ load x[6] from s6
eor r2, r2, r1    @ x[2] = x[2] ^ x[6]
vmov s2, r2    @ spill x[2] from r2
vmov r2, s28    @ load x[28] from s28
eor r10, r10, r2    @ x[13] = x[13] ^ x[28]
vmov s19, r0    @ spill x[19] from r0
vmov r0, s29    @ load x[29] from s29
eor r1, r1, r0    @ x[6] = x[6] ^ x[29]
eor r3, r3, r0    @ x[31] = x[31] ^ x[29]
vmov s29, r0    @ spill x[29] from r0
vmov r0, s5    @ load x[5] from s5
eor r4, r4, r0    @ x[1] = x[1] ^ x[5]
eor r0, r0, r3    @ x[5] = x[5] ^ x[31]
eor r5, r5, r8    @ x[7] = x[7] ^ x[30]
vmov s5, r0    @ spill x[5] from r0
vmov r0, s25    @ load x[25] from s25
eor r7, r7, r0    @ x[3] = x[3] ^ x[25]
vmov s3, r7    @ spill x[3] from r7
vmov r7, s0    @ load x[0] from s0
eor r7, r7, r0    @ x[0] = x[0] ^ x[25]
vmov s7, r5    @ spill x[7] from r5
vmov r5, s10    @ load x[10] from s10
vmov s30, r8    @ spill x[30] from r8
vmov r8, s20    @ load x[20] from s20
eor r5, r5, r8    @ x[10] = x[10] ^ x[20]
vmov s10, r5    @ spill x[10] from r5
vmov r5, s9    @ load x[9] from s9
eor r5, r5, r1    @ x[9] = x[9] ^ x[6]
eor r3, r3, r10    @ x[31] = x[31] ^ x[13]
eor r3, r3, r0    @ x[31] = x[31] ^ x[25]
eor r7, r7, r4    @ x[0] = x[0] ^ x[1]
eor r9, r9, r11    @ x[18] = x[18] ^ x[14]
vmov s31, r3    @ spill x[31] from r3
vmov r3, s22    @ load x[22] from s22
eor r3, r3, r2    @ x[22] = x[22] ^ x[28]
eor r2, r2, r6    @ x[28] = x[28] ^ x[26]
vmov s22, r3    @ spill x[22] from r3
vmov r3, s4    @ load x[4] from s4
eor r7, r7, r3    @ x[0] = x[0] ^ x[4]
eor r3, r3, r9    @ x[4] = x[4] ^ x[18]
vmov s18, r9    @ spill x[18] from r9
vmov r9, s19    @ load x[19] from s19
eor r0, r0, r9    @ x[25] = x[25] ^ x[19]
vmov s25, r0    @ spill x[25] from r0
vmov r0, s11    @ load x[11] from s11
eor r9, r9, r0    @ x[19] = x[19] ^ x[11]
vmov s19, r9    @ spill x[19] from r9
vmov r9, s24    @ load x[24] from s24
eor r1, r1, r9    @ x[6] = x[6] ^ x[24]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s29    @ load x[29] from s29
eor r1, r1, r0    @ x[29] = x[29] ^ x[11]
vmov s24, r9    @ spill x[24] from r9
vmov r9, s27    @ load x[27] from s27
eor r10, r10, r9    @ x[13] = x[13] ^ x[27]
eor r8, r8, r3    @ x[20] = x[20] ^ x[4]
vmov s27, r9    @ spill x[27] from r9
vmov r9, s30    @ load x[30] from s30
eor r4, r4, r9    @ x[1] = x[1] ^ x[30]
vmov s1, r4    @ spill x[1] from r4
vmov r4, s10    @ load x[10] from s10
eor r4, r4, r2    @ x[10] = x[10] ^ x[28]
vmov s20, r8    @ spill x[20] from r8
vmov r8, s17    @ load x[17] from s17
eor r8, r8, r7    @ x[17] = x[17] ^ x[0]
eor r3, r3, r9    @ x[4] = x[4] ^ x[30]
eor r11, r11, r9    @ x[14] = x[14] ^ x[30]
vmov s14, r11    @ spill x[14] from r11
vmov r11, s2    @ load x[2] from s2
eor r9, r9, r11    @ x[30] = x[30] ^ x[2]
vmov s2, r11    @ spill x[2] from r11
vmov r11, s21    @ load x[21] from s21
eor r11, r11, r6    @ x[21] = x[21] ^ x[26]
vmov s13, r10    @ spill x[13] from r10
vmov r10, s8    @ load x[8] from s8
eor r10, r10, r2    @ x[8] = x[8] ^ x[28]
vmov s28, r2    @ spill x[28] from r2
vmov r2, s12    @ load x[12] from s12
eor r2, r2, r1    @ x[12] = x[12] ^ x[29]
eor r1, r1, r4    @ x[29] = x[29] ^ x[10]
eor r5, r5, r7    @ x[9] = x[9] ^ x[0]
eor r11, r11, r3    @ x[21] = x[21] ^ x[4]
vmov s21, r11    @ spill x[21] from r11
vmov r11, s25    @ load x[25] from s25
eor r10, r10, r11    @ x[8] = x[8] ^ x[25]
vmov s9, r5    @ spill x[9] from r5
vmov r5, s31    @ load x[31] from s31
eor r5, r5, r6    @ x[31] = x[31] ^ x[26]
eor r4, r4, r8    @ x[10] = x[10] ^ x[17]
eor r11, r11, r9    @ x[25] = x[25] ^ x[30]
vmov s8, r10    @ spill x[8] from r10
vmov r10, s19    @ load x[19] from s19
eor r8, r8, r10    @ x[17] = x[17] ^ x[19]
vmov s31, r5    @ spill x[31] from r5
vmov r5, s18    @ load x[18] from s18
eor r5, r5, r2    @ x[18] = x[18] ^ x[12]
eor r9, r9, r0    @ x[30] = x[30] ^ x[11]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s24    @ load x[24] from s24
eor r0, r0, r9    @ x[11] = x[11] ^ x[24]
vmov s18, r5    @ spill x[18] from r5
vmov r5, s7    @ load x[7] from s7
eor r9, r9, r5    @ x[24] = x[24] ^ x[7]
vmov s11, r0    @ spill x[11] from r0
vmov r0, s3    @ load x[3] from s3
eor r9, r9, r0    @ x[24] = x[24] ^ x[3]
vmov s24, r9    @ spill x[24] from r9
vmov r9, s13    @ load x[13] from s13
eor r0, r0, r9    @ x[3] = x[3] ^ x[13]
eor r2, r2, r6    @ x[12] = x[12] ^ x[26]
vmov s25, r11    @ spill x[25] from r11
vmov r11, s20    @ load x[20] from s20
eor r6, r6, r11    @ x[26] = x[26] ^ x[20]
vmov s17, r8    @ spill x[17] from r8
vmov r8, s15    @ load x[15] from s15
eor r11, r11, r8    @ x[20] = x[20] ^ x[15]
eor r8, r8, r9    @ x[15] = x[15] ^ x[13]
eor r9, r9, r1    @ x[13] = x[13] ^ x[29]
eor r1, r1, r2    @ x[29] = x[29] ^ x[12]
eor r1, r1, r6    @ x[29] = x[29] ^ x[26]
eor r6, r6, r5    @ x[26] = x[26] ^ x[7]
vmov s12, r2    @ spill x[12] from r2
vmov r2, s6    @ load x[6] from s6
eor r5, r5, r2    @ x[7] = x[7] ^ x[6]
eor r8, r8, r2    @ x[15] = x[15] ^ x[6]
eor r11, r11, r10    @ x[20] = x[20] ^ x[19]
vmov s13, r9    @ spill x[13] from r9
vmov r9, s28    @ load x[28] from s28
eor r0, r0, r9    @ x[3] = x[3] ^ x[28]
vmov s15, r8    @ spill x[15] from r8
vmov r8, s2    @ load x[2] from s2
eor r2, r2, r8    @ x[6] = x[6] ^ x[2]
eor r10, r10, r9    @ x[19] = x[19] ^ x[28]
vmov s6, r2    @ spill x[6] from r2
vmov r2, s23    @ load x[23] from s23
eor r1, r1, r2    @ x[29] = x[29] ^ x[23]
vmov s29, r1    @ spill x[29] from r1
vmov r1, s22    @ load x[22] from s22
eor r5, r5, r1    @ x[7] = x[7] ^ x[22]
vmov s7, r5    @ spill x[7] from r5
vmov r5, s1    @ load x[1] from s1
eor r9, r9, r5    @ x[28] = x[28] ^ x[1]
vmov s23, r2    @ spill x[23] from r2
vmov r2, s5    @ load x[5] from s5
vmov s1, r5    @ spill x[1] from r5
vmov r5, s14    @ load x[14] from s14
eor r2, r2, r5    @ x[5] = x[5] ^ x[14]
eor r5, r5, r3    @ x[14] = x[14] ^ x[4]
eor r5, r5, r1    @ x[14] = x[14] ^ x[22]
vmov s5, r2    @ spill x[5] from r2
vmov r2, s16    @ load x[16] from s16
eor r9, r9, r2    @ x[28] = x[28] ^ x[16]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s27    @ load x[27] from s27
eor r11, r11, r9    @ x[20] = x[20] ^ x[27]
eor r10, r10, r3    @ x[19] = x[19] ^ x[4]
eor r7, r7, r1    @ x[0] = x[0] ^ x[22]
eor r3, r3, r4    @ x[4] = x[4] ^ x[10]
eor r4, r4, r8    @ x[10] = x[10] ^ x[2]
eor r8, r8, r0    @ x[2] = x[2] ^ x[3]
eor r8, r8, r6    @ x[2] = x[2] ^ x[26]
vmov s19, r10    @ spill x[19] from r10
vmov r10, s17    @ load x[17] from s17
vmov s10, r4    @ spill x[10] from r4
vmov r4, s15    @ load x[15] from s15
eor r10, r10, r4    @ x[17] = x[17] ^ x[15]
vmov s17, r10    @ spill x[17] from r10
vmov r10, s25    @ load x[25] from s25
eor r9, r9, r10    @ x[27] = x[27] ^ x[25]
eor r4, r4, r7    @ x[15] = x[15] ^ x[0]
eor r6, r6, r3    @ x[26] = x[26] ^ x[4]
vmov s0, r7    @ spill x[0] from r7
vmov r7, s31    @ load x[31] from s31
eor r3, r3, r7    @ x[4] = x[4] ^ x[31]
vmov s27, r9    @ spill x[27] from r9
vmov r9, s6    @ load x[6] from s6
vmov s31, r7    @ spill x[31] from r7
vmov r7, s11    @ load x[11] from s11
eor r9, r9, r7    @ x[6] = x[6] ^ x[11]
vmov s26, r6    @ spill x[26] from r6
vmov r6, s1    @ load x[1] from s1
vmov s15, r4    @ spill x[15] from r4
vmov r4, s23    @ load x[23] from s23
eor r6, r6, r4    @ x[1] = x[1] ^ x[23]
eor r4, r4, r9    @ x[23] = x[23] ^ x[6]
eor r10, r10, r8    @ x[25] = x[25] ^ x[2]
eor r9, r9, r3    @ x[6] = x[6] ^ x[4]
eor r10, r10, r1    @ x[25] = x[25] ^ x[22]
vmov s25, r10    @ spill x[25] from r10
vmov r10, s13    @ load x[13] from s13
eor r9, r9, r10    @ x[6] = x[6] ^ x[13]
vmov s6, r9    @ spill x[6] from r9
vmov r9, s18    @ load x[18] from s18
vmov s23, r4    @ spill x[23] from r4
vmov r4, s8    @ load x[8] from s8
eor r9, r9, r4    @ x[18] = x[18] ^ x[8]
eor r0, r0, r7    @ x[3] = x[3] ^ x[11]
eor r7, r7, r5    @ x[11] = x[11] ^ x[14]
vmov s18, r9    @ spill x[18] from r9
vmov r9, s9    @ load x[9] from s9
eor r0, r0, r9    @ x[3] = x[3] ^ x[9]
eor r10, r10, r2    @ x[13] = x[13] ^ x[16]
eor r2, r2, r8    @ x[16] = x[16] ^ x[2]
vmov s3, r0    @ spill x[3] from r0
vmov r0, s30    @ load x[30] from s30
eor r8, r8, r0    @ x[2] = x[2] ^ x[30]
eor r11, r11, r4    @ x[20] = x[20] ^ x[8]
vmov s20, r11    @ spill x[20] from r11
vmov r11, s10    @ load x[10] from s10
eor r11, r11, r4    @ x[10] = x[10] ^ x[8]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s29    @ load x[29] from s29
eor r1, r1, r11    @ x[22] = x[22] ^ x[29]
vmov s22, r1    @ spill x[22] from r1
vmov r1, s15    @ load x[15] from s15
eor r3, r3, r1    @ x[4] = x[4] ^ x[15]
vmov s16, r2    @ spill x[16] from r2
vmov r2, s21    @ load x[21] from s21
eor r4, r4, r2    @ x[8] = x[8] ^ x[21]
vmov s13, r10    @ spill x[13] from r10
vmov r10, s7    @ load x[7] from s7
eor r1, r1, r10    @ x[15] = x[15] ^ x[7]
eor r1, r1, r0    @ x[15] = x[15] ^ x[30]
eor r10, r10, r6    @ x[7] = x[7] ^ x[1]
vmov s7, r10    @ spill x[7] from r10
vmov r10, s26    @ load x[26] from s26
eor r6, r6, r10    @ x[1] = x[1] ^ x[26]
vmov s1, r6    @ spill x[1] from r6
vmov r6, s19    @ load x[19] from s19
eor r10, r10, r6    @ x[26] = x[26] ^ x[19]
eor r11, r11, r2    @ x[29] = x[29] ^ x[21]
eor r6, r6, r5    @ x[19] = x[19] ^ x[14]
eor r5, r5, r8    @ x[14] = x[14] ^ x[2]
vmov s29, r11    @ spill x[29] from r11
vmov r11, s5    @ load x[5] from s5
eor r0, r0, r11    @ x[30] = x[30] ^ x[5]
vmov s30, r0    @ spill x[30] from r0
vmov r0, s28    @ load x[28] from s28
eor r11, r11, r0    @ x[5] = x[5] ^ x[28]
vmov s28, r0    @ spill x[28] from r0
vmov r0, s23    @ load x[23] from s23
eor r8, r8, r0    @ x[2] = x[2] ^ x[23]
vmov s14, r5    @ spill x[14] from r5
vmov r5, s31    @ load x[31] from s31
eor r0, r0, r5    @ x[23] = x[23] ^ x[31]
eor r0, r0, r2    @ x[23] = x[23] ^ x[21]
vmov s23, r0    @ spill x[23] from r0
vmov r0, s12    @ load x[12] from s12
eor r6, r6, r0    @ x[19] = x[19] ^ x[12]
vmov s31, r5    @ spill x[31] from r5
vmov r5, s17    @ load x[17] from s17
eor r10, r10, r5    @ x[26] = x[26] ^ x[17]
eor r2, r2, r7    @ x[21] = x[21] ^ x[11]
eor r7, r7, r9    @ x[11] = x[11] ^ x[9]
eor r7, r7, r5    @ x[11] = x[11] ^ x[17]
eor r4, r4, r0    @ x[8] = x[8] ^ x[12]
eor r8, r8, r1    @ x[2] = x[2] ^ x[15]
vmov s15, r1    @ spill x[15] from r1
vmov r1, s6    @ load x[6] from s6
vmov s8, r4    @ spill x[8] from r4
vmov r4, s27    @ load x[27] from s27
eor r1, r1, r4    @ x[6] = x[6] ^ x[27]
eor r3, r3, r11    @ x[4] = x[4] ^ x[5]
eor r11, r11, r6    @ x[5] = x[5] ^ x[19]
eor r11, r11, r2    @ x[5] = x[5] ^ x[21]
eor r3, r3, r9    @ x[4] = x[4] ^ x[9]
eor r8, r8, r1    @ x[2] = x[2] ^ x[6]
eor r1, r1, r7    @ x[6] = x[6] ^ x[11]
vmov s5, r11    @ spill x[5] from r11
vmov r11, s13    @ load x[13] from s13
eor r5, r5, r11    @ x[17] = x[17] ^ x[13]
eor r2, r2, r11    @ x[21] = x[21] ^ x[13]
vmov s21, r2    @ spill x[21] from r2
vmov r2, s18    @ load x[18] from s18
eor r11, r11, r2    @ x[13] = x[13] ^ x[18]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s14    @ load x[14] from s14
eor r9, r9, r11    @ x[9] = x[9] ^ x[14]
eor r6, r6, r1    @ x[19] = x[19] ^ x[6]
vmov s2, r8    @ spill x[2] from r8
vmov r8, s31    @ load x[31] from s31
eor r0, r0, r8    @ x[12] = x[12] ^ x[31]
vmov s11, r7    @ spill x[11] from r7
vmov r7, s16    @ load x[16] from s16
eor r7, r7, r10    @ x[16] = x[16] ^ x[26]
eor r9, r9, r7    @ x[9] = x[9] ^ x[16]
eor r9, r9, r4    @ x[9] = x[9] ^ x[27]
vmov s9, r9    @ spill x[9] from r9
vmov r9, s30    @ load x[30] from s30
eor r9, r9, r8    @ x[30] = x[30] ^ x[31]
eor r8, r8, r2    @ x[31] = x[31] ^ x[18]
eor r7, r7, r5    @ x[16] = x[16] ^ x[17]
vmov s16, r7    @ spill x[16] from r7
vmov r7, s10    @ load x[10] from s10
eor r5, r5, r7    @ x[17] = x[17] ^ x[10]
vmov s17, r5    @ spill x[17] from r5
vmov r5, s0    @ load x[0] from s0
eor r1, r1, r5    @ x[6] = x[6] ^ x[0]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s8    @ load x[8] from s8
eor r1, r1, r5    @ x[8] = x[8] ^ x[0]
vmov s0, r5    @ spill x[0] from r5
vmov r5, s23    @ load x[23] from s23
eor r11, r11, r5    @ x[14] = x[14] ^ x[23]
eor r5, r5, r9    @ x[23] = x[23] ^ x[30]
eor r2, r2, r4    @ x[18] = x[18] ^ x[27]
eor r7, r7, r0    @ x[10] = x[10] ^ x[12]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s3    @ load x[3] from s3
eor r4, r4, r0    @ x[3] = x[3] ^ x[12]
eor r6, r6, r8    @ x[19] = x[19] ^ x[31]
eor r3, r3, r1    @ x[4] = x[4] ^ x[8]
eor r0, r0, r10    @ x[12] = x[12] ^ x[26]
eor r2, r2, r9    @ x[18] = x[18] ^ x[30]
vmov s26, r10    @ spill x[26] from r10
vmov r10, s17    @ load x[17] from s17
eor r10, r10, r9    @ x[17] = x[17] ^ x[30]
vmov s12, r0    @ spill x[12] from r0
vmov r0, s11    @ load x[11] from s11
eor r0, r0, r6    @ x[11] = x[11] ^ x[19]
vmov s18, r2    @ spill x[18] from r2
vmov r2, s15    @ load x[15] from s15
eor r6, r6, r2    @ x[19] = x[19] ^ x[15]
vmov s14, r11    @ spill x[14] from r11
vmov r11, s25    @ load x[25] from s25
eor r11, r11, r4    @ x[25] = x[25] ^ x[3]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s2    @ load x[2] from s2
eor r8, r8, r9    @ x[31] = x[31] ^ x[2]
vmov s10, r7    @ spill x[10] from r7
vmov r7, s20    @ load x[20] from s20
eor r7, r7, r3    @ x[20] = x[20] ^ x[4]
eor r1, r1, r9    @ x[8] = x[8] ^ x[2]
vmov s4, r3    @ spill x[4] from r3
vmov r3, s28    @ load x[28] from s28
vmov s8, r1    @ spill x[8] from r1
vmov r1, s16    @ load x[16] from s16
eor r3, r3, r1    @ x[28] = x[28] ^ x[16]
eor r11, r11, r3    @ x[25] = x[25] ^ x[28]
vmov s23, r5    @ spill x[23] from r5
vmov r5, s29    @ load x[29] from s29
eor r5, r5, r9    @ x[29] = x[29] ^ x[2]
vmov s3, r4    @ spill x[3] from r4
vmov r4, s1    @ load x[1] from s1
eor r4, r4, r6    @ x[1] = x[1] ^ x[19]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s22    @ load x[22] from s22
eor r6, r6, r7    @ x[22] = x[22] ^ x[20]
eor r8, r8, r6    @ x[31] = x[31] ^ x[22]
eor r10, r10, r4    @ x[17] = x[17] ^ x[1]
vstr.32  s6, [r14, #0]    @ y[0] = x[6] from s6
vstr.32  s19, [r14, #4]    @ y[1] = x[19] from s19
str  r2, [r14, #8]    @ y[2] = x[15] from r2
vstr.32  s12, [r14, #12]    @ y[3] = x[12] from s12
vstr.32  s18, [r14, #16]    @ y[4] = x[18] from s18
vstr.32  s3, [r14, #20]    @ y[5] = x[3] from s3
vstr.32  s26, [r14, #24]    @ y[6] = x[26] from s26
vstr.32  s30, [r14, #28]    @ y[7] = x[30] from s30
vstr.32  s10, [r14, #32]    @ y[8] = x[10] from s10
str  r7, [r14, #36]    @ y[9] = x[20] from r7
vstr.32  s9, [r14, #40]    @ y[10] = x[9] from s9
str  r4, [r14, #44]    @ y[11] = x[1] from r4
str  r0, [r14, #48]    @ y[12] = x[11] from r0
str  r10, [r14, #52]    @ y[13] = x[17] from r10
vstr.32  s8, [r14, #56]    @ y[14] = x[8] from s8
str  r11, [r14, #60]    @ y[15] = x[25] from r11
vstr.32  s27, [r14, #64]    @ y[16] = x[27] from s27
vstr.32  s24, [r14, #68]    @ y[17] = x[24] from s24
str  r1, [r14, #72]    @ y[18] = x[16] from r1
str  r3, [r14, #76]    @ y[19] = x[28] from r3
vstr.32  s7, [r14, #80]    @ y[20] = x[7] from s7
vstr.32  s0, [r14, #84]    @ y[21] = x[0] from s0
vstr.32  s23, [r14, #88]    @ y[22] = x[23] from s23
vstr.32  s14, [r14, #92]    @ y[23] = x[14] from s14
vstr.32  s21, [r14, #96]    @ y[24] = x[21] from s21
str  r9, [r14, #100]    @ y[25] = x[2] from r9
vstr.32  s4, [r14, #104]    @ y[26] = x[4] from s4
vstr.32  s13, [r14, #108]    @ y[27] = x[13] from s13
str  r6, [r14, #112]    @ y[28] = x[22] from r6
str  r8, [r14, #116]    @ y[29] = x[31] from r8
vstr.32  s5, [r14, #120]    @ y[30] = x[5] from s5
str  r5, [r14, #124]    @ y[31] = x[29] from r5
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v21, .-gft_mul_v21
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v22
.type gft_mul_v22, %function
.align 2
gft_mul_v22:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s2    @ load x[2] from s2
vmov r1, s3    @ load x[3] from s3
eor r0, r0, r1    @ x[2] = x[2] ^ x[3]
vmov r2, s29    @ load x[29] from s29
eor r0, r0, r2    @ x[2] = x[2] ^ x[29]
vmov r3, s4    @ load x[4] from s4
eor r3, r3, r0    @ x[4] = x[4] ^ x[2]
vmov r4, s24    @ load x[24] from s24
eor r4, r4, r3    @ x[24] = x[24] ^ x[4]
vmov r5, s10    @ load x[10] from s10
vmov r6, s9    @ load x[9] from s9
eor r5, r5, r6    @ x[10] = x[10] ^ x[9]
vmov r7, s18    @ load x[18] from s18
vmov r8, s23    @ load x[23] from s23
eor r7, r7, r8    @ x[18] = x[18] ^ x[23]
vmov r9, s11    @ load x[11] from s11
eor r9, r9, r5    @ x[11] = x[11] ^ x[10]
vmov r10, s26    @ load x[26] from s26
eor r10, r10, r4    @ x[26] = x[26] ^ x[24]
vmov r11, s28    @ load x[28] from s28
eor r11, r11, r0    @ x[28] = x[28] ^ x[2]
vmov s9, r6    @ spill x[9] from r6
vmov r6, s8    @ load x[8] from s8
eor r6, r6, r9    @ x[8] = x[8] ^ x[11]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s20    @ load x[20] from s20
vmov s26, r10    @ spill x[26] from r10
vmov r10, s22    @ load x[22] from s22
eor r4, r4, r10    @ x[20] = x[20] ^ x[22]
eor r10, r10, r7    @ x[22] = x[22] ^ x[18]
vmov s18, r7    @ spill x[18] from r7
vmov r7, s1    @ load x[1] from s1
eor r3, r3, r7    @ x[4] = x[4] ^ x[1]
vmov s22, r10    @ spill x[22] from r10
vmov r10, s27    @ load x[27] from s27
vmov s2, r0    @ spill x[2] from r0
vmov r0, s6    @ load x[6] from s6
eor r10, r10, r0    @ x[27] = x[27] ^ x[6]
vmov s27, r10    @ spill x[27] from r10
vmov r10, s5    @ load x[5] from s5
vmov s6, r0    @ spill x[6] from r0
vmov r0, s12    @ load x[12] from s12
eor r10, r10, r0    @ x[5] = x[5] ^ x[12]
vmov s5, r10    @ spill x[5] from r10
vmov r10, s16    @ load x[16] from s16
eor r10, r10, r6    @ x[16] = x[16] ^ x[8]
eor r6, r6, r7    @ x[8] = x[8] ^ x[1]
vmov s16, r10    @ spill x[16] from r10
vmov r10, s17    @ load x[17] from s17
eor r9, r9, r10    @ x[11] = x[11] ^ x[17]
vmov s11, r9    @ spill x[11] from r9
vmov r9, s21    @ load x[21] from s21
vmov s8, r6    @ spill x[8] from r6
vmov r6, s31    @ load x[31] from s31
eor r9, r9, r6    @ x[21] = x[21] ^ x[31]
eor r1, r1, r3    @ x[3] = x[3] ^ x[4]
vmov s31, r6    @ spill x[31] from r6
vmov r6, s0    @ load x[0] from s0
eor r5, r5, r6    @ x[10] = x[10] ^ x[0]
vmov s10, r5    @ spill x[10] from r5
vmov r5, s7    @ load x[7] from s7
eor r10, r10, r5    @ x[17] = x[17] ^ x[7]
eor r8, r8, r5    @ x[23] = x[23] ^ x[7]
eor r10, r10, r11    @ x[17] = x[17] ^ x[28]
eor r11, r11, r4    @ x[28] = x[28] ^ x[20]
eor r6, r6, r4    @ x[0] = x[0] ^ x[20]
eor r4, r4, r2    @ x[20] = x[20] ^ x[29]
vmov s28, r11    @ spill x[28] from r11
vmov r11, s2    @ load x[2] from s2
vmov s23, r8    @ spill x[23] from r8
vmov r8, s13    @ load x[13] from s13
eor r11, r11, r8    @ x[2] = x[2] ^ x[13]
vmov s20, r4    @ spill x[20] from r4
vmov r4, s25    @ load x[25] from s25
eor r2, r2, r4    @ x[29] = x[29] ^ x[25]
eor r8, r8, r2    @ x[13] = x[13] ^ x[29]
vmov s13, r8    @ spill x[13] from r8
vmov r8, s26    @ load x[26] from s26
eor r5, r5, r8    @ x[7] = x[7] ^ x[26]
vmov s2, r11    @ spill x[2] from r11
vmov r11, s30    @ load x[30] from s30
eor r8, r8, r11    @ x[26] = x[26] ^ x[30]
eor r2, r2, r7    @ x[29] = x[29] ^ x[1]
vmov s26, r8    @ spill x[26] from r8
vmov r8, s15    @ load x[15] from s15
vmov s25, r4    @ spill x[25] from r4
vmov r4, s14    @ load x[14] from s14
eor r8, r8, r4    @ x[15] = x[15] ^ x[14]
eor r1, r1, r0    @ x[3] = x[3] ^ x[12]
vmov s15, r8    @ spill x[15] from r8
vmov r8, s22    @ load x[22] from s22
eor r11, r11, r8    @ x[30] = x[30] ^ x[22]
vmov s14, r4    @ spill x[14] from r4
vmov r4, s6    @ load x[6] from s6
vmov s3, r1    @ spill x[3] from r1
vmov r1, s24    @ load x[24] from s24
eor r4, r4, r1    @ x[6] = x[6] ^ x[24]
vmov s6, r4    @ spill x[6] from r4
vmov r4, s8    @ load x[8] from s8
eor r4, r4, r5    @ x[8] = x[8] ^ x[7]
vmov s7, r5    @ spill x[7] from r5
vmov r5, s5    @ load x[5] from s5
eor r5, r5, r4    @ x[5] = x[5] ^ x[8]
eor r4, r4, r10    @ x[8] = x[8] ^ x[17]
vmov s8, r4    @ spill x[8] from r4
vmov r4, s27    @ load x[27] from s27
eor r5, r5, r4    @ x[5] = x[5] ^ x[27]
eor r10, r10, r11    @ x[17] = x[17] ^ x[30]
eor r11, r11, r9    @ x[30] = x[30] ^ x[21]
vmov s17, r10    @ spill x[17] from r10
vmov r10, s19    @ load x[19] from s19
eor r7, r7, r10    @ x[1] = x[1] ^ x[19]
eor r2, r2, r0    @ x[29] = x[29] ^ x[12]
eor r9, r9, r1    @ x[21] = x[21] ^ x[24]
eor r10, r10, r6    @ x[19] = x[19] ^ x[0]
eor r0, r0, r1    @ x[12] = x[12] ^ x[24]
eor r6, r6, r3    @ x[0] = x[0] ^ x[4]
vmov s29, r2    @ spill x[29] from r2
vmov r2, s10    @ load x[10] from s10
eor r2, r2, r1    @ x[10] = x[10] ^ x[24]
eor r3, r3, r8    @ x[4] = x[4] ^ x[22]
vmov s0, r6    @ spill x[0] from r6
vmov r6, s25    @ load x[25] from s25
eor r4, r4, r6    @ x[27] = x[27] ^ x[25]
vmov s21, r9    @ spill x[21] from r9
vmov r9, s31    @ load x[31] from s31
eor r8, r8, r9    @ x[22] = x[22] ^ x[31]
vmov s4, r3    @ spill x[4] from r3
vmov r3, s2    @ load x[2] from s2
vmov s22, r8    @ spill x[22] from r8
vmov r8, s6    @ load x[6] from s6
eor r3, r3, r8    @ x[2] = x[2] ^ x[6]
vmov s2, r3    @ spill x[2] from r3
vmov r3, s20    @ load x[20] from s20
eor r3, r3, r0    @ x[20] = x[20] ^ x[12]
vmov s12, r0    @ spill x[12] from r0
vmov r0, s16    @ load x[16] from s16
eor r0, r0, r6    @ x[16] = x[16] ^ x[25]
vmov s16, r0    @ spill x[16] from r0
vmov r0, s11    @ load x[11] from s11
vmov s5, r5    @ spill x[5] from r5
vmov r5, s7    @ load x[7] from s7
eor r0, r0, r5    @ x[11] = x[11] ^ x[7]
eor r5, r5, r10    @ x[7] = x[7] ^ x[19]
eor r10, r10, r2    @ x[19] = x[19] ^ x[10]
vmov s19, r10    @ spill x[19] from r10
vmov r10, s18    @ load x[18] from s18
eor r6, r6, r10    @ x[25] = x[25] ^ x[18]
eor r1, r1, r5    @ x[24] = x[24] ^ x[7]
vmov s7, r5    @ spill x[7] from r5
vmov r5, s8    @ load x[8] from s8
eor r10, r10, r5    @ x[18] = x[18] ^ x[8]
eor r5, r5, r2    @ x[8] = x[8] ^ x[10]
eor r2, r2, r4    @ x[10] = x[10] ^ x[27]
eor r5, r5, r9    @ x[8] = x[8] ^ x[31]
vmov s8, r5    @ spill x[8] from r5
vmov r5, s3    @ load x[3] from s3
eor r5, r5, r8    @ x[3] = x[3] ^ x[6]
eor r1, r1, r6    @ x[24] = x[24] ^ x[25]
eor r6, r6, r4    @ x[25] = x[25] ^ x[27]
vmov s24, r1    @ spill x[24] from r1
vmov r1, s14    @ load x[14] from s14
eor r1, r1, r9    @ x[14] = x[14] ^ x[31]
eor r6, r6, r0    @ x[25] = x[25] ^ x[11]
eor r9, r9, r4    @ x[31] = x[31] ^ x[27]
eor r0, r0, r10    @ x[11] = x[11] ^ x[18]
vmov s25, r6    @ spill x[25] from r6
vmov r6, s23    @ load x[23] from s23
eor r9, r9, r6    @ x[31] = x[31] ^ x[23]
eor r6, r6, r7    @ x[23] = x[23] ^ x[1]
vmov s10, r2    @ spill x[10] from r2
vmov r2, s15    @ load x[15] from s15
eor r0, r0, r2    @ x[11] = x[11] ^ x[15]
vmov s31, r9    @ spill x[31] from r9
vmov r9, s9    @ load x[9] from s9
eor r7, r7, r9    @ x[1] = x[1] ^ x[9]
eor r10, r10, r11    @ x[18] = x[18] ^ x[30]
eor r9, r9, r10    @ x[9] = x[9] ^ x[18]
eor r4, r4, r2    @ x[27] = x[27] ^ x[15]
eor r10, r10, r3    @ x[18] = x[18] ^ x[20]
vmov s1, r7    @ spill x[1] from r7
vmov r7, s5    @ load x[5] from s5
eor r2, r2, r7    @ x[15] = x[15] ^ x[5]
eor r3, r3, r7    @ x[20] = x[20] ^ x[5]
vmov s20, r3    @ spill x[20] from r3
vmov r3, s19    @ load x[19] from s19
eor r10, r10, r3    @ x[18] = x[18] ^ x[19]
vmov s14, r1    @ spill x[14] from r1
vmov r1, s28    @ load x[28] from s28
eor r7, r7, r1    @ x[5] = x[5] ^ x[28]
vmov s11, r0    @ spill x[11] from r0
vmov r0, s16    @ load x[16] from s16
eor r8, r8, r0    @ x[6] = x[6] ^ x[16]
eor r7, r7, r3    @ x[5] = x[5] ^ x[19]
vmov s6, r8    @ spill x[6] from r8
vmov r8, s2    @ load x[2] from s2
eor r2, r2, r8    @ x[15] = x[15] ^ x[2]
eor r3, r3, r4    @ x[19] = x[19] ^ x[27]
vmov s19, r3    @ spill x[19] from r3
vmov r3, s22    @ load x[22] from s22
eor r4, r4, r3    @ x[27] = x[27] ^ x[22]
eor r6, r6, r1    @ x[23] = x[23] ^ x[28]
eor r10, r10, r8    @ x[18] = x[18] ^ x[2]
eor r3, r3, r1    @ x[22] = x[22] ^ x[28]
vmov s18, r10    @ spill x[18] from r10
vmov r10, s7    @ load x[7] from s7
eor r1, r1, r10    @ x[28] = x[28] ^ x[7]
eor r9, r9, r4    @ x[9] = x[9] ^ x[27]
eor r8, r8, r4    @ x[2] = x[2] ^ x[27]
vmov s9, r9    @ spill x[9] from r9
vmov r9, s4    @ load x[4] from s4
eor r4, r4, r9    @ x[27] = x[27] ^ x[4]
eor r7, r7, r5    @ x[5] = x[5] ^ x[3]
eor r10, r10, r11    @ x[7] = x[7] ^ x[30]
eor r11, r11, r0    @ x[30] = x[30] ^ x[16]
eor r4, r4, r5    @ x[27] = x[27] ^ x[3]
vmov s5, r7    @ spill x[5] from r7
vmov r7, s13    @ load x[13] from s13
eor r3, r3, r7    @ x[22] = x[22] ^ x[13]
eor r5, r5, r7    @ x[3] = x[3] ^ x[13]
eor r0, r0, r9    @ x[16] = x[16] ^ x[4]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s21    @ load x[21] from s21
eor r9, r9, r4    @ x[4] = x[4] ^ x[21]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s11    @ load x[11] from s11
eor r6, r6, r11    @ x[23] = x[23] ^ x[11]
eor r4, r4, r11    @ x[21] = x[21] ^ x[11]
vmov s21, r4    @ spill x[21] from r4
vmov r4, s26    @ load x[26] from s26
eor r7, r7, r4    @ x[13] = x[13] ^ x[26]
vmov s16, r0    @ spill x[16] from r0
vmov r0, s0    @ load x[0] from s0
eor r3, r3, r0    @ x[22] = x[22] ^ x[0]
vmov s22, r3    @ spill x[22] from r3
vmov r3, s6    @ load x[6] from s6
eor r11, r11, r3    @ x[11] = x[11] ^ x[6]
eor r10, r10, r0    @ x[7] = x[7] ^ x[0]
eor r2, r2, r7    @ x[15] = x[15] ^ x[13]
vmov s15, r2    @ spill x[15] from r2
vmov r2, s14    @ load x[14] from s14
eor r7, r7, r2    @ x[13] = x[13] ^ x[14]
eor r4, r4, r1    @ x[26] = x[26] ^ x[28]
eor r4, r4, r10    @ x[26] = x[26] ^ x[7]
eor r6, r6, r2    @ x[23] = x[23] ^ x[14]
eor r10, r10, r5    @ x[7] = x[7] ^ x[3]
eor r5, r5, r2    @ x[3] = x[3] ^ x[14]
vmov s23, r6    @ spill x[23] from r6
vmov r6, s31    @ load x[31] from s31
eor r8, r8, r6    @ x[2] = x[2] ^ x[31]
eor r9, r9, r6    @ x[4] = x[4] ^ x[31]
eor r6, r6, r10    @ x[31] = x[31] ^ x[7]
eor r2, r2, r3    @ x[14] = x[14] ^ x[6]
vmov s14, r2    @ spill x[14] from r2
vmov r2, s19    @ load x[19] from s19
eor r10, r10, r2    @ x[7] = x[7] ^ x[19]
eor r2, r2, r3    @ x[19] = x[19] ^ x[6]
eor r3, r3, r0    @ x[6] = x[6] ^ x[0]
vmov s13, r7    @ spill x[13] from r7
vmov r7, s12    @ load x[12] from s12
eor r0, r0, r7    @ x[0] = x[0] ^ x[12]
vmov s6, r3    @ spill x[6] from r3
vmov r3, s29    @ load x[29] from s29
eor r11, r11, r3    @ x[11] = x[11] ^ x[29]
vmov s31, r6    @ spill x[31] from r6
vmov r6, s16    @ load x[16] from s16
eor r7, r7, r6    @ x[12] = x[12] ^ x[16]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s17    @ load x[17] from s17
eor r1, r1, r2    @ x[28] = x[28] ^ x[17]
vmov s28, r1    @ spill x[28] from r1
vmov r1, s20    @ load x[20] from s20
eor r1, r1, r6    @ x[20] = x[20] ^ x[16]
eor r6, r6, r2    @ x[16] = x[16] ^ x[17]
vmov s3, r5    @ spill x[3] from r5
vmov r5, s30    @ load x[30] from s30
vmov s7, r10    @ spill x[7] from r10
vmov r10, s10    @ load x[10] from s10
eor r5, r5, r10    @ x[30] = x[30] ^ x[10]
vmov s10, r10    @ spill x[10] from r10
vmov r10, s1    @ load x[1] from s1
eor r7, r7, r10    @ x[12] = x[12] ^ x[1]
vmov s12, r7    @ spill x[12] from r7
vmov r7, s25    @ load x[25] from s25
eor r6, r6, r7    @ x[16] = x[16] ^ x[25]
vmov s25, r7    @ spill x[25] from r7
vmov r7, s27    @ load x[27] from s27
eor r7, r7, r6    @ x[27] = x[27] ^ x[16]
eor r1, r1, r9    @ x[20] = x[20] ^ x[4]
eor r3, r3, r10    @ x[29] = x[29] ^ x[1]
vmov s20, r1    @ spill x[20] from r1
vmov r1, s9    @ load x[9] from s9
vmov s29, r3    @ spill x[29] from r3
vmov r3, s8    @ load x[8] from s8
eor r1, r1, r3    @ x[9] = x[9] ^ x[8]
eor r3, r3, r2    @ x[8] = x[8] ^ x[17]
vmov s8, r3    @ spill x[8] from r3
vmov r3, s22    @ load x[22] from s22
eor r3, r3, r5    @ x[22] = x[22] ^ x[30]
eor r8, r8, r11    @ x[2] = x[2] ^ x[11]
eor r11, r11, r7    @ x[11] = x[11] ^ x[27]
eor r7, r7, r0    @ x[27] = x[27] ^ x[0]
eor r0, r0, r4    @ x[0] = x[0] ^ x[26]
vmov s2, r8    @ spill x[2] from r8
vmov r8, s7    @ load x[7] from s7
eor r8, r8, r0    @ x[7] = x[7] ^ x[0]
vmov s27, r7    @ spill x[27] from r7
vmov r7, s5    @ load x[5] from s5
eor r10, r10, r7    @ x[1] = x[1] ^ x[5]
vmov s7, r8    @ spill x[7] from r8
vmov r8, s3    @ load x[3] from s3
eor r8, r8, r10    @ x[3] = x[3] ^ x[1]
vmov s1, r10    @ spill x[1] from r10
vmov r10, s18    @ load x[18] from s18
eor r10, r10, r8    @ x[18] = x[18] ^ x[3]
vmov s3, r8    @ spill x[3] from r8
vmov r8, s28    @ load x[28] from s28
eor r4, r4, r8    @ x[26] = x[26] ^ x[28]
eor r6, r6, r4    @ x[16] = x[16] ^ x[26]
eor r2, r2, r9    @ x[17] = x[17] ^ x[4]
vmov s18, r10    @ spill x[18] from r10
vmov r10, s25    @ load x[25] from s25
eor r10, r10, r8    @ x[25] = x[25] ^ x[28]
vmov s28, r8    @ spill x[28] from r8
vmov r8, s21    @ load x[21] from s21
eor r8, r8, r1    @ x[21] = x[21] ^ x[9]
vmov s16, r6    @ spill x[16] from r6
vmov r6, s10    @ load x[10] from s10
eor r1, r1, r6    @ x[9] = x[9] ^ x[10]
eor r5, r5, r9    @ x[30] = x[30] ^ x[4]
eor r5, r5, r7    @ x[30] = x[30] ^ x[5]
vmov s4, r9    @ spill x[4] from r9
vmov r9, s19    @ load x[19] from s19
eor r7, r7, r9    @ x[5] = x[5] ^ x[19]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s31    @ load x[31] from s31
vmov s19, r9    @ spill x[19] from r9
vmov r9, s8    @ load x[8] from s8
eor r6, r6, r9    @ x[31] = x[31] ^ x[8]
eor r7, r7, r1    @ x[5] = x[5] ^ x[9]
vmov s9, r1    @ spill x[9] from r1
vmov r1, s6    @ load x[6] from s6
eor r1, r1, r11    @ x[6] = x[6] ^ x[11]
vmov s11, r11    @ spill x[11] from r11
vmov r11, s13    @ load x[13] from s13
eor r11, r11, r8    @ x[13] = x[13] ^ x[21]
eor r3, r3, r2    @ x[22] = x[22] ^ x[17]
eor r0, r0, r10    @ x[0] = x[0] ^ x[25]
eor r4, r4, r5    @ x[26] = x[26] ^ x[30]
vmov s26, r4    @ spill x[26] from r4
vmov r4, s14    @ load x[14] from s14
vmov s0, r0    @ spill x[0] from r0
vmov r0, s16    @ load x[16] from s16
eor r4, r4, r0    @ x[14] = x[14] ^ x[16]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s29    @ load x[29] from s29
eor r2, r2, r3    @ x[29] = x[29] ^ x[22]
vstr.32  s15, [r14, #0]    @ y[0] = x[15] from s15
str  r11, [r14, #4]    @ y[1] = x[13] from r11
vstr.32  s7, [r14, #8]    @ y[2] = x[7] from s7
vstr.32  s3, [r14, #12]    @ y[3] = x[3] from s3
vstr.32  s11, [r14, #16]    @ y[4] = x[11] from s11
str  r8, [r14, #20]    @ y[5] = x[21] from r8
vstr.32  s10, [r14, #24]    @ y[6] = x[10] from s10
vstr.32  s1, [r14, #28]    @ y[7] = x[1] from s1
str  r6, [r14, #32]    @ y[8] = x[31] from r6
vstr.32  s0, [r14, #36]    @ y[9] = x[0] from s0
vstr.32  s27, [r14, #40]    @ y[10] = x[27] from s27
vstr.32  s18, [r14, #44]    @ y[11] = x[18] from s18
str  r9, [r14, #48]    @ y[12] = x[8] from r9
str  r10, [r14, #52]    @ y[13] = x[25] from r10
str  r0, [r14, #56]    @ y[14] = x[16] from r0
str  r7, [r14, #60]    @ y[15] = x[5] from r7
str  r5, [r14, #64]    @ y[16] = x[30] from r5
vstr.32  s28, [r14, #68]    @ y[17] = x[28] from s28
str  r1, [r14, #72]    @ y[18] = x[6] from r1
vstr.32  s17, [r14, #76]    @ y[19] = x[17] from s17
vstr.32  s26, [r14, #80]    @ y[20] = x[26] from s26
vstr.32  s20, [r14, #84]    @ y[21] = x[20] from s20
vstr.32  s24, [r14, #88]    @ y[22] = x[24] from s24
vstr.32  s4, [r14, #92]    @ y[23] = x[4] from s4
str  r2, [r14, #96]    @ y[24] = x[29] from r2
vstr.32  s12, [r14, #100]    @ y[25] = x[12] from s12
vstr.32  s19, [r14, #104]    @ y[26] = x[19] from s19
str  r4, [r14, #108]    @ y[27] = x[14] from r4
str  r3, [r14, #112]    @ y[28] = x[22] from r3
vstr.32  s2, [r14, #116]    @ y[29] = x[2] from s2
vstr.32  s23, [r14, #120]    @ y[30] = x[23] from s23
vstr.32  s9, [r14, #124]    @ y[31] = x[9] from s9
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v22, .-gft_mul_v22
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v23
.type gft_mul_v23, %function
.align 2
gft_mul_v23:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s21    @ load x[21] from s21
vmov r1, s18    @ load x[18] from s18
eor r0, r0, r1    @ x[21] = x[21] ^ x[18]
vmov r2, s2    @ load x[2] from s2
vmov r3, s24    @ load x[24] from s24
eor r2, r2, r3    @ x[2] = x[2] ^ x[24]
vmov r4, s16    @ load x[16] from s16
eor r4, r4, r0    @ x[16] = x[16] ^ x[21]
vmov r5, s23    @ load x[23] from s23
eor r5, r5, r4    @ x[23] = x[23] ^ x[16]
vmov r6, s22    @ load x[22] from s22
eor r6, r6, r5    @ x[22] = x[22] ^ x[23]
vmov r7, s19    @ load x[19] from s19
eor r6, r6, r7    @ x[22] = x[22] ^ x[19]
vmov r8, s5    @ load x[5] from s5
eor r8, r8, r6    @ x[5] = x[5] ^ x[22]
vmov r9, s26    @ load x[26] from s26
eor r1, r1, r9    @ x[18] = x[18] ^ x[26]
eor r9, r9, r2    @ x[26] = x[26] ^ x[2]
vmov r10, s8    @ load x[8] from s8
vmov r11, s30    @ load x[30] from s30
eor r10, r10, r11    @ x[8] = x[8] ^ x[30]
eor r0, r0, r2    @ x[21] = x[21] ^ x[2]
eor r3, r3, r0    @ x[24] = x[24] ^ x[21]
vmov s23, r5    @ spill x[23] from r5
vmov r5, s14    @ load x[14] from s14
vmov s22, r6    @ spill x[22] from r6
vmov r6, s25    @ load x[25] from s25
eor r5, r5, r6    @ x[14] = x[14] ^ x[25]
vmov s16, r4    @ spill x[16] from r4
vmov r4, s4    @ load x[4] from s4
eor r4, r4, r7    @ x[4] = x[4] ^ x[19]
vmov s19, r7    @ spill x[19] from r7
vmov r7, s1    @ load x[1] from s1
eor r7, r7, r9    @ x[1] = x[1] ^ x[26]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s3    @ load x[3] from s3
eor r11, r11, r8    @ x[3] = x[3] ^ x[5]
eor r7, r7, r10    @ x[1] = x[1] ^ x[8]
vmov s1, r7    @ spill x[1] from r7
vmov r7, s12    @ load x[12] from s12
eor r9, r9, r7    @ x[26] = x[26] ^ x[12]
vmov s3, r11    @ spill x[3] from r11
vmov r11, s11    @ load x[11] from s11
eor r11, r11, r9    @ x[11] = x[11] ^ x[26]
vmov s11, r11    @ spill x[11] from r11
vmov r11, s17    @ load x[17] from s17
eor r0, r0, r11    @ x[21] = x[21] ^ x[17]
vmov s21, r0    @ spill x[21] from r0
vmov r0, s0    @ load x[0] from s0
eor r0, r0, r10    @ x[0] = x[0] ^ x[8]
vmov s17, r11    @ spill x[17] from r11
vmov r11, s31    @ load x[31] from s31
eor r1, r1, r11    @ x[18] = x[18] ^ x[31]
eor r6, r6, r8    @ x[25] = x[25] ^ x[5]
eor r10, r10, r3    @ x[8] = x[8] ^ x[24]
vmov s8, r10    @ spill x[8] from r10
vmov r10, s6    @ load x[6] from s6
eor r10, r10, r4    @ x[6] = x[6] ^ x[4]
vmov s25, r6    @ spill x[25] from r6
vmov r6, s9    @ load x[9] from s9
eor r6, r6, r5    @ x[9] = x[9] ^ x[14]
vmov s9, r6    @ spill x[9] from r6
vmov r6, s27    @ load x[27] from s27
eor r6, r6, r2    @ x[27] = x[27] ^ x[2]
vmov s31, r11    @ spill x[31] from r11
vmov r11, s28    @ load x[28] from s28
eor r11, r11, r7    @ x[28] = x[28] ^ x[12]
vmov s6, r10    @ spill x[6] from r10
vmov r10, s30    @ load x[30] from s30
eor r7, r7, r10    @ x[12] = x[12] ^ x[30]
vmov s12, r7    @ spill x[12] from r7
vmov r7, s19    @ load x[19] from s19
eor r10, r10, r7    @ x[30] = x[30] ^ x[19]
vmov s30, r10    @ spill x[30] from r10
vmov r10, s3    @ load x[3] from s3
eor r7, r7, r10    @ x[19] = x[19] ^ x[3]
eor r10, r10, r9    @ x[3] = x[3] ^ x[26]
vmov s27, r6    @ spill x[27] from r6
vmov r6, s16    @ load x[16] from s16
eor r7, r7, r6    @ x[19] = x[19] ^ x[16]
vmov s18, r1    @ spill x[18] from r1
vmov r1, s29    @ load x[29] from s29
eor r2, r2, r1    @ x[2] = x[2] ^ x[29]
eor r0, r0, r2    @ x[0] = x[0] ^ x[2]
eor r6, r6, r2    @ x[16] = x[16] ^ x[2]
eor r7, r7, r3    @ x[19] = x[19] ^ x[24]
vmov s19, r7    @ spill x[19] from r7
vmov r7, s22    @ load x[22] from s22
eor r3, r3, r7    @ x[24] = x[24] ^ x[22]
eor r2, r2, r7    @ x[2] = x[2] ^ x[22]
eor r6, r6, r4    @ x[16] = x[16] ^ x[4]
vmov s16, r6    @ spill x[16] from r6
vmov r6, s10    @ load x[10] from s10
eor r7, r7, r6    @ x[22] = x[22] ^ x[10]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s11    @ load x[11] from s11
eor r6, r6, r7    @ x[10] = x[10] ^ x[11]
eor r7, r7, r1    @ x[11] = x[11] ^ x[29]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s7    @ load x[7] from s7
vmov s24, r3    @ spill x[24] from r3
vmov r3, s1    @ load x[1] from s1
eor r6, r6, r3    @ x[7] = x[7] ^ x[1]
eor r4, r4, r8    @ x[4] = x[4] ^ x[5]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s17    @ load x[17] from s17
eor r6, r6, r11    @ x[17] = x[17] ^ x[28]
eor r10, r10, r6    @ x[3] = x[3] ^ x[17]
eor r3, r3, r1    @ x[1] = x[1] ^ x[29]
eor r9, r9, r5    @ x[26] = x[26] ^ x[14]
eor r5, r5, r8    @ x[14] = x[14] ^ x[5]
vmov s28, r11    @ spill x[28] from r11
vmov r11, s23    @ load x[23] from s23
eor r6, r6, r11    @ x[17] = x[17] ^ x[23]
eor r6, r6, r1    @ x[17] = x[17] ^ x[29]
vmov s14, r5    @ spill x[14] from r5
vmov r5, s20    @ load x[20] from s20
eor r2, r2, r5    @ x[2] = x[2] ^ x[20]
eor r6, r6, r2    @ x[17] = x[17] ^ x[2]
vmov s17, r6    @ spill x[17] from r6
vmov r6, s18    @ load x[18] from s18
eor r7, r7, r6    @ x[11] = x[11] ^ x[18]
vmov s3, r10    @ spill x[3] from r10
vmov r10, s12    @ load x[12] from s12
vmov s29, r1    @ spill x[29] from r1
vmov r1, s27    @ load x[27] from s27
eor r10, r10, r1    @ x[12] = x[12] ^ x[27]
eor r0, r0, r8    @ x[0] = x[0] ^ x[5]
vmov s26, r9    @ spill x[26] from r9
vmov r9, s24    @ load x[24] from s24
eor r8, r8, r9    @ x[5] = x[5] ^ x[24]
eor r7, r7, r8    @ x[11] = x[11] ^ x[5]
vmov s11, r7    @ spill x[11] from r7
vmov r7, s30    @ load x[30] from s30
eor r8, r8, r7    @ x[5] = x[5] ^ x[30]
eor r7, r7, r10    @ x[30] = x[30] ^ x[12]
eor r5, r5, r6    @ x[20] = x[20] ^ x[18]
eor r6, r6, r3    @ x[18] = x[18] ^ x[1]
eor r6, r6, r4    @ x[18] = x[18] ^ x[4]
vmov s20, r5    @ spill x[20] from r5
vmov r5, s7    @ load x[7] from s7
eor r0, r0, r5    @ x[0] = x[0] ^ x[7]
vmov s30, r7    @ spill x[30] from r7
vmov r7, s19    @ load x[19] from s19
eor r7, r7, r3    @ x[19] = x[19] ^ x[1]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s10    @ load x[10] from s10
eor r8, r8, r6    @ x[5] = x[5] ^ x[10]
vmov s19, r7    @ spill x[19] from r7
vmov r7, s6    @ load x[6] from s6
eor r7, r7, r0    @ x[6] = x[6] ^ x[0]
eor r11, r11, r6    @ x[23] = x[23] ^ x[10]
vmov s6, r7    @ spill x[6] from r7
vmov r7, s22    @ load x[22] from s22
eor r4, r4, r7    @ x[4] = x[4] ^ x[22]
vmov s1, r3    @ spill x[1] from r3
vmov r3, s13    @ load x[13] from s13
eor r6, r6, r3    @ x[10] = x[10] ^ x[13]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s21    @ load x[21] from s21
eor r2, r2, r6    @ x[2] = x[2] ^ x[21]
eor r7, r7, r1    @ x[22] = x[22] ^ x[27]
eor r1, r1, r6    @ x[27] = x[27] ^ x[21]
eor r4, r4, r9    @ x[4] = x[4] ^ x[24]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s31    @ load x[31] from s31
eor r11, r11, r7    @ x[23] = x[23] ^ x[31]
vmov s4, r4    @ spill x[4] from r4
vmov r4, s25    @ load x[25] from s25
eor r3, r3, r4    @ x[13] = x[13] ^ x[25]
eor r9, r9, r0    @ x[24] = x[24] ^ x[0]
eor r0, r0, r4    @ x[0] = x[0] ^ x[25]
eor r4, r4, r5    @ x[25] = x[25] ^ x[7]
eor r6, r6, r5    @ x[21] = x[21] ^ x[7]
vmov s25, r4    @ spill x[25] from r4
vmov r4, s26    @ load x[26] from s26
eor r5, r5, r4    @ x[7] = x[7] ^ x[26]
eor r4, r4, r8    @ x[26] = x[26] ^ x[5]
eor r2, r2, r7    @ x[2] = x[2] ^ x[31]
vmov s2, r2    @ spill x[2] from r2
vmov r2, s29    @ load x[29] from s29
eor r2, r2, r10    @ x[29] = x[29] ^ x[12]
vmov s29, r2    @ spill x[29] from r2
vmov r2, s1    @ load x[1] from s1
vmov s26, r4    @ spill x[26] from r4
vmov r4, s15    @ load x[15] from s15
eor r2, r2, r4    @ x[1] = x[1] ^ x[15]
vmov s1, r2    @ spill x[1] from r2
vmov r2, s3    @ load x[3] from s3
eor r6, r6, r2    @ x[21] = x[21] ^ x[3]
vmov s21, r6    @ spill x[21] from r6
vmov r6, s14    @ load x[14] from s14
eor r6, r6, r1    @ x[14] = x[14] ^ x[27]
vmov s24, r9    @ spill x[24] from r9
vmov r9, s28    @ load x[28] from s28
eor r1, r1, r9    @ x[27] = x[27] ^ x[28]
eor r1, r1, r0    @ x[27] = x[27] ^ x[0]
eor r8, r8, r0    @ x[5] = x[5] ^ x[0]
eor r0, r0, r3    @ x[0] = x[0] ^ x[13]
vmov s0, r0    @ spill x[0] from r0
vmov r0, s8    @ load x[8] from s8
eor r3, r3, r0    @ x[13] = x[13] ^ x[8]
eor r10, r10, r7    @ x[12] = x[12] ^ x[31]
eor r3, r3, r11    @ x[13] = x[13] ^ x[23]
eor r6, r6, r4    @ x[14] = x[14] ^ x[15]
eor r4, r4, r0    @ x[15] = x[15] ^ x[8]
eor r1, r1, r2    @ x[27] = x[27] ^ x[3]
eor r5, r5, r7    @ x[7] = x[7] ^ x[31]
vmov s3, r2    @ spill x[3] from r2
vmov r2, s24    @ load x[24] from s24
eor r4, r4, r2    @ x[15] = x[15] ^ x[24]
vmov s15, r4    @ spill x[15] from r4
vmov r4, s17    @ load x[17] from s17
eor r4, r4, r1    @ x[17] = x[17] ^ x[27]
eor r7, r7, r0    @ x[31] = x[31] ^ x[8]
eor r7, r7, r9    @ x[31] = x[31] ^ x[28]
vmov s27, r1    @ spill x[27] from r1
vmov r1, s4    @ load x[4] from s4
eor r8, r8, r1    @ x[5] = x[5] ^ x[4]
vmov s31, r7    @ spill x[31] from r7
vmov r7, s16    @ load x[16] from s16
eor r9, r9, r7    @ x[28] = x[28] ^ x[16]
vmov s5, r8    @ spill x[5] from r8
vmov r8, s21    @ load x[21] from s21
vmov s4, r1    @ spill x[4] from r1
vmov r1, s19    @ load x[19] from s19
eor r8, r8, r1    @ x[21] = x[21] ^ x[19]
vmov s21, r8    @ spill x[21] from r8
vmov r8, s26    @ load x[26] from s26
eor r8, r8, r1    @ x[26] = x[26] ^ x[19]
vmov s26, r8    @ spill x[26] from r8
vmov r8, s18    @ load x[18] from s18
eor r1, r1, r8    @ x[19] = x[19] ^ x[18]
vmov s18, r8    @ spill x[18] from r8
vmov r8, s30    @ load x[30] from s30
vmov s16, r7    @ spill x[16] from r7
vmov r7, s20    @ load x[20] from s20
eor r8, r8, r7    @ x[30] = x[30] ^ x[20]
vmov s30, r8    @ spill x[30] from r8
vmov r8, s11    @ load x[11] from s11
eor r8, r8, r3    @ x[11] = x[11] ^ x[13]
eor r2, r2, r4    @ x[24] = x[24] ^ x[17]
eor r0, r0, r5    @ x[8] = x[8] ^ x[7]
vmov s11, r8    @ spill x[11] from r8
vmov r8, s9    @ load x[9] from s9
eor r5, r5, r8    @ x[7] = x[7] ^ x[9]
vmov s7, r5    @ spill x[7] from r5
vmov r5, s6    @ load x[6] from s6
eor r9, r9, r5    @ x[28] = x[28] ^ x[6]
eor r11, r11, r1    @ x[23] = x[23] ^ x[19]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s22    @ load x[22] from s22
eor r9, r9, r10    @ x[22] = x[22] ^ x[12]
eor r3, r3, r6    @ x[13] = x[13] ^ x[14]
eor r6, r6, r5    @ x[14] = x[14] ^ x[6]
vmov s13, r3    @ spill x[13] from r3
vmov r3, s2    @ load x[2] from s2
eor r0, r0, r3    @ x[8] = x[8] ^ x[2]
eor r1, r1, r4    @ x[19] = x[19] ^ x[17]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s16    @ load x[16] from s16
eor r2, r2, r6    @ x[24] = x[24] ^ x[16]
eor r7, r7, r8    @ x[20] = x[20] ^ x[9]
vmov s19, r1    @ spill x[19] from r1
vmov r1, s0    @ load x[0] from s0
eor r8, r8, r1    @ x[9] = x[9] ^ x[0]
eor r1, r1, r10    @ x[0] = x[0] ^ x[12]
vmov s0, r1    @ spill x[0] from r1
vmov r1, s25    @ load x[25] from s25
eor r10, r10, r1    @ x[12] = x[12] ^ x[25]
vmov s24, r2    @ spill x[24] from r2
vmov r2, s1    @ load x[1] from s1
eor r1, r1, r2    @ x[25] = x[25] ^ x[1]
vmov s20, r7    @ spill x[20] from r7
vmov r7, s7    @ load x[7] from s7
eor r2, r2, r7    @ x[1] = x[1] ^ x[7]
vmov s22, r9    @ spill x[22] from r9
vmov r9, s29    @ load x[29] from s29
eor r7, r7, r9    @ x[7] = x[7] ^ x[29]
eor r9, r9, r11    @ x[29] = x[29] ^ x[23]
vmov s29, r9    @ spill x[29] from r9
vmov r9, s4    @ load x[4] from s4
eor r11, r11, r9    @ x[23] = x[23] ^ x[4]
eor r9, r9, r5    @ x[4] = x[4] ^ x[6]
eor r5, r5, r2    @ x[6] = x[6] ^ x[1]
vmov s4, r9    @ spill x[4] from r9
vmov r9, s5    @ load x[5] from s5
eor r1, r1, r9    @ x[25] = x[25] ^ x[5]
vmov s25, r1    @ spill x[25] from r1
vmov r1, s3    @ load x[3] from s3
eor r6, r6, r1    @ x[16] = x[16] ^ x[3]
vmov s16, r6    @ spill x[16] from r6
vmov r6, s15    @ load x[15] from s15
eor r8, r8, r6    @ x[9] = x[9] ^ x[15]
eor r10, r10, r1    @ x[12] = x[12] ^ x[3]
eor r9, r9, r8    @ x[5] = x[5] ^ x[9]
eor r8, r8, r3    @ x[9] = x[9] ^ x[2]
vmov s9, r8    @ spill x[9] from r8
vmov r8, s30    @ load x[30] from s30
eor r3, r3, r8    @ x[2] = x[2] ^ x[30]
vmov s5, r9    @ spill x[5] from r9
vmov r9, s10    @ load x[10] from s10
eor r3, r3, r9    @ x[2] = x[2] ^ x[10]
eor r4, r4, r7    @ x[17] = x[17] ^ x[7]
eor r5, r5, r4    @ x[6] = x[6] ^ x[17]
vmov s17, r4    @ spill x[17] from r4
vmov r4, s28    @ load x[28] from s28
eor r0, r0, r4    @ x[8] = x[8] ^ x[28]
vmov s8, r0    @ spill x[8] from r0
vmov r0, s26    @ load x[26] from s26
eor r0, r0, r4    @ x[26] = x[26] ^ x[28]
eor r11, r11, r1    @ x[23] = x[23] ^ x[3]
vmov s23, r11    @ spill x[23] from r11
vmov r11, s31    @ load x[31] from s31
eor r10, r10, r11    @ x[12] = x[12] ^ x[31]
eor r7, r7, r9    @ x[7] = x[7] ^ x[10]
eor r4, r4, r8    @ x[28] = x[28] ^ x[30]
eor r5, r5, r6    @ x[6] = x[6] ^ x[15]
eor r9, r9, r0    @ x[10] = x[10] ^ x[26]
vmov s3, r1    @ spill x[3] from r1
vmov r1, s21    @ load x[21] from s21
eor r1, r1, r2    @ x[21] = x[21] ^ x[1]
vmov s28, r4    @ spill x[28] from r4
vmov r4, s4    @ load x[4] from s4
eor r11, r11, r4    @ x[31] = x[31] ^ x[4]
eor r9, r9, r11    @ x[10] = x[10] ^ x[31]
vmov s10, r9    @ spill x[10] from r9
vmov r9, s22    @ load x[22] from s22
vmov s26, r0    @ spill x[26] from r0
vmov r0, s20    @ load x[20] from s20
eor r9, r9, r0    @ x[22] = x[22] ^ x[20]
vmov s31, r11    @ spill x[31] from r11
vmov r11, s24    @ load x[24] from s24
eor r10, r10, r11    @ x[12] = x[12] ^ x[24]
eor r0, r0, r3    @ x[20] = x[20] ^ x[2]
vmov s2, r3    @ spill x[2] from r3
vmov r3, s11    @ load x[11] from s11
eor r3, r3, r1    @ x[11] = x[11] ^ x[21]
eor r6, r6, r7    @ x[15] = x[15] ^ x[7]
vmov s7, r7    @ spill x[7] from r7
vmov r7, s29    @ load x[29] from s29
eor r7, r7, r10    @ x[29] = x[29] ^ x[12]
vmov s12, r10    @ spill x[12] from r10
vmov r10, s27    @ load x[27] from s27
eor r10, r10, r5    @ x[27] = x[27] ^ x[6]
eor r8, r8, r5    @ x[30] = x[30] ^ x[6]
vmov s6, r5    @ spill x[6] from r5
vmov r5, s16    @ load x[16] from s16
eor r5, r5, r10    @ x[16] = x[16] ^ x[27]
vmov s15, r6    @ spill x[15] from r6
vmov r6, s5    @ load x[5] from s5
eor r2, r2, r6    @ x[1] = x[1] ^ x[5]
vstr.32  s10, [r14, #0]    @ y[0] = x[10] from s10
vstr.32  s25, [r14, #4]    @ y[1] = x[25] from s25
vstr.32  s7, [r14, #8]    @ y[2] = x[7] from s7
vstr.32  s31, [r14, #12]    @ y[3] = x[31] from s31
str  r6, [r14, #16]    @ y[4] = x[5] from r6
vstr.32  s26, [r14, #20]    @ y[5] = x[26] from s26
str  r2, [r14, #24]    @ y[6] = x[1] from r2
vstr.32  s15, [r14, #28]    @ y[7] = x[15] from s15
vstr.32  s17, [r14, #32]    @ y[8] = x[17] from s17
vstr.32  s18, [r14, #36]    @ y[9] = x[18] from s18
str  r7, [r14, #40]    @ y[10] = x[29] from r7
str  r4, [r14, #44]    @ y[11] = x[4] from r4
vstr.32  s0, [r14, #48]    @ y[12] = x[0] from s0
vstr.32  s6, [r14, #52]    @ y[13] = x[6] from s6
vstr.32  s14, [r14, #56]    @ y[14] = x[14] from s14
vstr.32  s8, [r14, #60]    @ y[15] = x[8] from s8
vstr.32  s12, [r14, #64]    @ y[16] = x[12] from s12
vstr.32  s28, [r14, #68]    @ y[17] = x[28] from s28
str  r11, [r14, #72]    @ y[18] = x[24] from r11
vstr.32  s3, [r14, #76]    @ y[19] = x[3] from s3
str  r8, [r14, #80]    @ y[20] = x[30] from r8
str  r10, [r14, #84]    @ y[21] = x[27] from r10
str  r5, [r14, #88]    @ y[22] = x[16] from r5
vstr.32  s19, [r14, #92]    @ y[23] = x[19] from s19
str  r1, [r14, #96]    @ y[24] = x[21] from r1
vstr.32  s9, [r14, #100]    @ y[25] = x[9] from s9
vstr.32  s23, [r14, #104]    @ y[26] = x[23] from s23
vstr.32  s2, [r14, #108]    @ y[27] = x[2] from s2
vstr.32  s13, [r14, #112]    @ y[28] = x[13] from s13
str  r9, [r14, #116]    @ y[29] = x[22] from r9
str  r0, [r14, #120]    @ y[30] = x[20] from r0
str  r3, [r14, #124]    @ y[31] = x[11] from r3
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v23, .-gft_mul_v23
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v24
.type gft_mul_v24, %function
.align 2
gft_mul_v24:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s22    @ load x[22] from s22
vmov r1, s26    @ load x[26] from s26
eor r0, r0, r1    @ x[22] = x[22] ^ x[26]
vmov r2, s30    @ load x[30] from s30
vmov r3, s31    @ load x[31] from s31
eor r2, r2, r3    @ x[30] = x[30] ^ x[31]
vmov r4, s6    @ load x[6] from s6
vmov r5, s15    @ load x[15] from s15
eor r4, r4, r5    @ x[6] = x[6] ^ x[15]
vmov r6, s12    @ load x[12] from s12
vmov r7, s13    @ load x[13] from s13
eor r6, r6, r7    @ x[12] = x[12] ^ x[13]
vmov r8, s2    @ load x[2] from s2
eor r8, r8, r4    @ x[2] = x[2] ^ x[6]
vmov r9, s23    @ load x[23] from s23
vmov r10, s10    @ load x[10] from s10
eor r9, r9, r10    @ x[23] = x[23] ^ x[10]
vmov r11, s11    @ load x[11] from s11
vmov s31, r3    @ spill x[31] from r3
vmov r3, s1    @ load x[1] from s1
eor r11, r11, r3    @ x[11] = x[11] ^ x[1]
vmov s13, r7    @ spill x[13] from r7
vmov r7, s5    @ load x[5] from s5
vmov s6, r4    @ spill x[6] from r4
vmov r4, s24    @ load x[24] from s24
eor r7, r7, r4    @ x[5] = x[5] ^ x[24]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s29    @ load x[29] from s29
eor r4, r4, r0    @ x[29] = x[29] ^ x[22]
eor r0, r0, r2    @ x[22] = x[22] ^ x[30]
vmov s30, r2    @ spill x[30] from r2
vmov r2, s27    @ load x[27] from s27
eor r2, r2, r9    @ x[27] = x[27] ^ x[23]
eor r1, r1, r6    @ x[26] = x[26] ^ x[12]
vmov s27, r2    @ spill x[27] from r2
vmov r2, s18    @ load x[18] from s18
eor r2, r2, r0    @ x[18] = x[18] ^ x[22]
vmov s22, r0    @ spill x[22] from r0
vmov r0, s16    @ load x[16] from s16
eor r5, r5, r0    @ x[15] = x[15] ^ x[16]
eor r6, r6, r11    @ x[12] = x[12] ^ x[11]
vmov s12, r6    @ spill x[12] from r6
vmov r6, s25    @ load x[25] from s25
vmov s15, r5    @ spill x[15] from r5
vmov r5, s21    @ load x[21] from s21
eor r6, r6, r5    @ x[25] = x[25] ^ x[21]
eor r10, r10, r0    @ x[10] = x[10] ^ x[16]
eor r5, r5, r7    @ x[21] = x[21] ^ x[5]
vmov s10, r10    @ spill x[10] from r10
vmov r10, s3    @ load x[3] from s3
eor r11, r11, r10    @ x[11] = x[11] ^ x[3]
eor r0, r0, r8    @ x[16] = x[16] ^ x[2]
vmov s11, r11    @ spill x[11] from r11
vmov r11, s20    @ load x[20] from s20
eor r7, r7, r11    @ x[5] = x[5] ^ x[20]
vmov s20, r11    @ spill x[20] from r11
vmov r11, s28    @ load x[28] from s28
eor r7, r7, r11    @ x[5] = x[5] ^ x[28]
eor r1, r1, r11    @ x[26] = x[26] ^ x[28]
vmov s28, r11    @ spill x[28] from r11
vmov r11, s9    @ load x[9] from s9
vmov s25, r6    @ spill x[25] from r6
vmov r6, s17    @ load x[17] from s17
eor r11, r11, r6    @ x[9] = x[9] ^ x[17]
vmov s16, r0    @ spill x[16] from r0
vmov r0, s19    @ load x[19] from s19
eor r9, r9, r0    @ x[23] = x[23] ^ x[19]
vmov s5, r7    @ spill x[5] from r7
vmov r7, s8    @ load x[8] from s8
eor r3, r3, r7    @ x[1] = x[1] ^ x[8]
vmov s8, r7    @ spill x[8] from r7
vmov r7, s6    @ load x[6] from s6
eor r7, r7, r2    @ x[6] = x[6] ^ x[18]
vmov s6, r7    @ spill x[6] from r7
vmov r7, s0    @ load x[0] from s0
eor r2, r2, r7    @ x[18] = x[18] ^ x[0]
vmov s18, r2    @ spill x[18] from r2
vmov r2, s14    @ load x[14] from s14
eor r4, r4, r2    @ x[29] = x[29] ^ x[14]
vmov s29, r4    @ spill x[29] from r4
vmov r4, s22    @ load x[22] from s22
eor r4, r4, r10    @ x[22] = x[22] ^ x[3]
vmov s0, r7    @ spill x[0] from r7
vmov r7, s4    @ load x[4] from s4
eor r10, r10, r7    @ x[3] = x[3] ^ x[4]
vmov s4, r7    @ spill x[4] from r7
vmov r7, s24    @ load x[24] from s24
eor r6, r6, r7    @ x[17] = x[17] ^ x[24]
eor r4, r4, r0    @ x[22] = x[22] ^ x[19]
eor r8, r8, r11    @ x[2] = x[2] ^ x[9]
eor r5, r5, r3    @ x[21] = x[21] ^ x[1]
eor r9, r9, r1    @ x[23] = x[23] ^ x[26]
vmov s22, r4    @ spill x[22] from r4
vmov r4, s5    @ load x[5] from s5
eor r8, r8, r4    @ x[2] = x[2] ^ x[5]
vmov s23, r9    @ spill x[23] from r9
vmov r9, s30    @ load x[30] from s30
eor r4, r4, r9    @ x[5] = x[5] ^ x[30]
eor r1, r1, r4    @ x[26] = x[26] ^ x[5]
eor r4, r4, r3    @ x[5] = x[5] ^ x[1]
eor r8, r8, r5    @ x[2] = x[2] ^ x[21]
vmov s26, r1    @ spill x[26] from r1
vmov r1, s7    @ load x[7] from s7
eor r5, r5, r1    @ x[21] = x[21] ^ x[7]
vmov s21, r5    @ spill x[21] from r5
vmov r5, s15    @ load x[15] from s15
eor r4, r4, r5    @ x[5] = x[5] ^ x[15]
eor r5, r5, r10    @ x[15] = x[15] ^ x[3]
vmov s5, r4    @ spill x[5] from r4
vmov r4, s8    @ load x[8] from s8
eor r5, r5, r4    @ x[15] = x[15] ^ x[8]
eor r8, r8, r0    @ x[2] = x[2] ^ x[19]
eor r10, r10, r11    @ x[3] = x[3] ^ x[9]
eor r10, r10, r9    @ x[3] = x[3] ^ x[30]
vmov s2, r8    @ spill x[2] from r8
vmov r8, s27    @ load x[27] from s27
eor r4, r4, r8    @ x[8] = x[8] ^ x[27]
vmov s8, r4    @ spill x[8] from r4
vmov r4, s16    @ load x[16] from s16
eor r9, r9, r4    @ x[30] = x[30] ^ x[16]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s25    @ load x[25] from s25
eor r4, r4, r9    @ x[16] = x[16] ^ x[25]
eor r9, r9, r11    @ x[25] = x[25] ^ x[9]
eor r11, r11, r8    @ x[9] = x[9] ^ x[27]
vmov s25, r9    @ spill x[25] from r9
vmov r9, s28    @ load x[28] from s28
eor r8, r8, r9    @ x[27] = x[27] ^ x[28]
vmov s27, r8    @ spill x[27] from r8
vmov r8, s20    @ load x[20] from s20
eor r4, r4, r8    @ x[16] = x[16] ^ x[20]
eor r11, r11, r2    @ x[9] = x[9] ^ x[14]
eor r2, r2, r0    @ x[14] = x[14] ^ x[19]
vmov s16, r4    @ spill x[16] from r4
vmov r4, s11    @ load x[11] from s11
eor r9, r9, r4    @ x[28] = x[28] ^ x[11]
eor r0, r0, r4    @ x[19] = x[19] ^ x[11]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s10    @ load x[10] from s10
eor r8, r8, r9    @ x[20] = x[20] ^ x[10]
eor r9, r9, r6    @ x[10] = x[10] ^ x[17]
eor r4, r4, r9    @ x[11] = x[11] ^ x[10]
vmov s11, r4    @ spill x[11] from r4
vmov r4, s0    @ load x[0] from s0
eor r2, r2, r4    @ x[14] = x[14] ^ x[0]
vmov s9, r11    @ spill x[9] from r11
vmov r11, s13    @ load x[13] from s13
eor r4, r4, r11    @ x[0] = x[0] ^ x[13]
vmov s0, r4    @ spill x[0] from r4
vmov r4, s18    @ load x[18] from s18
eor r9, r9, r4    @ x[10] = x[10] ^ x[18]
eor r4, r4, r7    @ x[18] = x[18] ^ x[24]
vmov s10, r9    @ spill x[10] from r9
vmov r9, s12    @ load x[12] from s12
eor r11, r11, r9    @ x[13] = x[13] ^ x[12]
eor r9, r9, r3    @ x[12] = x[12] ^ x[1]
eor r11, r11, r10    @ x[13] = x[13] ^ x[3]
eor r10, r10, r7    @ x[3] = x[3] ^ x[24]
eor r1, r1, r5    @ x[7] = x[7] ^ x[15]
eor r5, r5, r7    @ x[15] = x[15] ^ x[24]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s31    @ load x[31] from s31
eor r3, r3, r11    @ x[1] = x[1] ^ x[31]
eor r0, r0, r7    @ x[19] = x[19] ^ x[24]
eor r7, r7, r11    @ x[24] = x[24] ^ x[31]
eor r2, r2, r11    @ x[14] = x[14] ^ x[31]
eor r11, r11, r6    @ x[31] = x[31] ^ x[17]
vmov s19, r0    @ spill x[19] from r0
vmov r0, s29    @ load x[29] from s29
eor r5, r5, r0    @ x[15] = x[15] ^ x[29]
vmov s1, r3    @ spill x[1] from r3
vmov r3, s6    @ load x[6] from s6
eor r9, r9, r3    @ x[12] = x[12] ^ x[6]
vmov s7, r1    @ spill x[7] from r1
vmov r1, s4    @ load x[4] from s4
eor r6, r6, r1    @ x[17] = x[17] ^ x[4]
vmov s6, r3    @ spill x[6] from r3
vmov r3, s25    @ load x[25] from s25
eor r0, r0, r3    @ x[29] = x[29] ^ x[25]
vmov s3, r10    @ spill x[3] from r10
vmov r10, s27    @ load x[27] from s27
eor r8, r8, r10    @ x[20] = x[20] ^ x[27]
vmov s14, r2    @ spill x[14] from r2
vmov r2, s26    @ load x[26] from s26
eor r2, r2, r3    @ x[26] = x[26] ^ x[25]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s21    @ load x[21] from s21
eor r7, r7, r3    @ x[24] = x[24] ^ x[21]
eor r1, r1, r3    @ x[4] = x[4] ^ x[21]
eor r3, r3, r9    @ x[21] = x[21] ^ x[12]
vmov s21, r3    @ spill x[21] from r3
vmov r3, s8    @ load x[8] from s8
eor r9, r9, r3    @ x[12] = x[12] ^ x[8]
vmov s12, r9    @ spill x[12] from r9
vmov r9, s0    @ load x[0] from s0
eor r9, r9, r2    @ x[0] = x[0] ^ x[26]
vmov s8, r3    @ spill x[8] from r3
vmov r3, s9    @ load x[9] from s9
eor r9, r9, r3    @ x[0] = x[0] ^ x[9]
vmov s9, r3    @ spill x[9] from r3
vmov r3, s16    @ load x[16] from s16
eor r3, r3, r9    @ x[16] = x[16] ^ x[0]
vmov s0, r9    @ spill x[0] from r9
vmov r9, s23    @ load x[23] from s23
eor r9, r9, r11    @ x[23] = x[23] ^ x[31]
vmov s23, r9    @ spill x[23] from r9
vmov r9, s28    @ load x[28] from s28
eor r9, r9, r4    @ x[28] = x[28] ^ x[18]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s30    @ load x[30] from s30
eor r1, r1, r9    @ x[4] = x[4] ^ x[30]
vmov s4, r1    @ spill x[4] from r1
vmov r1, s11    @ load x[11] from s11
eor r1, r1, r5    @ x[11] = x[11] ^ x[15]
eor r10, r10, r3    @ x[27] = x[27] ^ x[16]
eor r6, r6, r8    @ x[17] = x[17] ^ x[20]
eor r4, r4, r2    @ x[18] = x[18] ^ x[26]
eor r2, r2, r8    @ x[26] = x[26] ^ x[20]
eor r8, r8, r7    @ x[20] = x[20] ^ x[24]
eor r6, r6, r0    @ x[17] = x[17] ^ x[29]
eor r8, r8, r3    @ x[20] = x[20] ^ x[16]
vmov s17, r6    @ spill x[17] from r6
vmov r6, s8    @ load x[8] from s8
eor r3, r3, r6    @ x[16] = x[16] ^ x[8]
eor r3, r3, r11    @ x[16] = x[16] ^ x[31]
eor r6, r6, r9    @ x[8] = x[8] ^ x[30]
vmov s18, r4    @ spill x[18] from r4
vmov r4, s25    @ load x[25] from s25
vmov s27, r10    @ spill x[27] from r10
vmov r10, s12    @ load x[12] from s12
eor r4, r4, r10    @ x[25] = x[25] ^ x[12]
vmov s26, r2    @ spill x[26] from r2
vmov r2, s14    @ load x[14] from s14
eor r0, r0, r2    @ x[29] = x[29] ^ x[14]
vmov s29, r0    @ spill x[29] from r0
vmov r0, s21    @ load x[21] from s21
vmov s25, r4    @ spill x[25] from r4
vmov r4, s10    @ load x[10] from s10
eor r0, r0, r4    @ x[21] = x[21] ^ x[10]
eor r4, r4, r2    @ x[10] = x[10] ^ x[14]
vmov s21, r0    @ spill x[21] from r0
vmov r0, s3    @ load x[3] from s3
eor r11, r11, r0    @ x[31] = x[31] ^ x[3]
vmov s11, r1    @ spill x[11] from r1
vmov r1, s6    @ load x[6] from s6
eor r9, r9, r1    @ x[30] = x[30] ^ x[6]
eor r2, r2, r9    @ x[14] = x[14] ^ x[30]
eor r9, r9, r5    @ x[30] = x[30] ^ x[15]
eor r9, r9, r7    @ x[30] = x[30] ^ x[24]
eor r4, r4, r0    @ x[10] = x[10] ^ x[3]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s7    @ load x[7] from s7
eor r0, r0, r1    @ x[3] = x[3] ^ x[7]
vmov s14, r2    @ spill x[14] from r2
vmov r2, s23    @ load x[23] from s23
eor r5, r5, r2    @ x[15] = x[15] ^ x[23]
vmov s15, r5    @ spill x[15] from r5
vmov r5, s1    @ load x[1] from s1
eor r6, r6, r5    @ x[8] = x[8] ^ x[1]
vmov s8, r6    @ spill x[8] from r6
vmov r6, s28    @ load x[28] from s28
eor r10, r10, r6    @ x[12] = x[12] ^ x[28]
vmov s1, r5    @ spill x[1] from r5
vmov r5, s9    @ load x[9] from s9
eor r11, r11, r5    @ x[31] = x[31] ^ x[9]
vmov s31, r11    @ spill x[31] from r11
vmov r11, s22    @ load x[22] from s22
eor r4, r4, r11    @ x[10] = x[10] ^ x[22]
vmov s10, r4    @ spill x[10] from r4
vmov r4, s0    @ load x[0] from s0
eor r10, r10, r4    @ x[12] = x[12] ^ x[0]
vmov s22, r11    @ spill x[22] from r11
vmov r11, s19    @ load x[19] from s19
eor r6, r6, r11    @ x[28] = x[28] ^ x[19]
vmov s19, r11    @ spill x[19] from r11
vmov r11, s2    @ load x[2] from s2
eor r7, r7, r11    @ x[24] = x[24] ^ x[2]
vmov s24, r7    @ spill x[24] from r7
vmov r7, s5    @ load x[5] from s5
eor r4, r4, r7    @ x[0] = x[0] ^ x[5]
eor r8, r8, r5    @ x[20] = x[20] ^ x[9]
eor r3, r3, r7    @ x[16] = x[16] ^ x[5]
eor r10, r10, r11    @ x[12] = x[12] ^ x[2]
eor r6, r6, r5    @ x[28] = x[28] ^ x[9]
vmov s16, r3    @ spill x[16] from r3
vmov r3, s13    @ load x[13] from s13
eor r0, r0, r3    @ x[3] = x[3] ^ x[13]
eor r7, r7, r9    @ x[5] = x[5] ^ x[30]
eor r9, r9, r1    @ x[30] = x[30] ^ x[7]
vmov s5, r7    @ spill x[5] from r7
vmov r7, s14    @ load x[14] from s14
eor r7, r7, r2    @ x[14] = x[14] ^ x[23]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s6    @ load x[6] from s6
vmov s23, r2    @ spill x[23] from r2
vmov r2, s19    @ load x[19] from s19
eor r9, r9, r2    @ x[6] = x[6] ^ x[19]
vmov s28, r6    @ spill x[28] from r6
vmov r6, s22    @ load x[22] from s22
eor r1, r1, r6    @ x[7] = x[7] ^ x[22]
vmov s3, r0    @ spill x[3] from r0
vmov r0, s4    @ load x[4] from s4
eor r11, r11, r0    @ x[2] = x[2] ^ x[4]
vmov s13, r3    @ spill x[13] from r3
vmov r3, s10    @ load x[10] from s10
eor r10, r10, r3    @ x[12] = x[12] ^ x[10]
vmov s22, r6    @ spill x[22] from r6
vmov r6, s11    @ load x[11] from s11
eor r3, r3, r6    @ x[10] = x[10] ^ x[11]
vmov s11, r6    @ spill x[11] from r6
vmov r6, s25    @ load x[25] from s25
eor r5, r5, r6    @ x[9] = x[9] ^ x[25]
vmov s10, r3    @ spill x[10] from r3
vmov r3, s26    @ load x[26] from s26
eor r6, r6, r3    @ x[25] = x[25] ^ x[26]
vmov s26, r3    @ spill x[26] from r3
vmov r3, s31    @ load x[31] from s31
eor r3, r3, r9    @ x[31] = x[31] ^ x[6]
eor r4, r4, r1    @ x[0] = x[0] ^ x[7]
vmov s31, r3    @ spill x[31] from r3
vmov r3, s21    @ load x[21] from s21
eor r3, r3, r5    @ x[21] = x[21] ^ x[9]
eor r1, r1, r11    @ x[7] = x[7] ^ x[2]
eor r2, r2, r8    @ x[19] = x[19] ^ x[20]
eor r0, r0, r3    @ x[4] = x[4] ^ x[21]
vmov s7, r1    @ spill x[7] from r1
vmov r1, s27    @ load x[27] from s27
eor r9, r9, r1    @ x[6] = x[6] ^ x[27]
eor r7, r7, r10    @ x[14] = x[14] ^ x[12]
vmov s6, r9    @ spill x[6] from r9
vmov r9, s1    @ load x[1] from s1
eor r2, r2, r9    @ x[19] = x[19] ^ x[1]
eor r11, r11, r8    @ x[2] = x[2] ^ x[20]
eor r8, r8, r6    @ x[20] = x[20] ^ x[25]
vmov s4, r0    @ spill x[4] from r0
vmov r0, s22    @ load x[22] from s22
eor r6, r6, r0    @ x[25] = x[25] ^ x[22]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s15    @ load x[15] from s15
eor r0, r0, r2    @ x[22] = x[22] ^ x[15]
vmov s27, r1    @ spill x[27] from r1
vmov r1, s13    @ load x[13] from s13
eor r6, r6, r1    @ x[25] = x[25] ^ x[13]
eor r2, r2, r9    @ x[15] = x[15] ^ x[1]
eor r1, r1, r9    @ x[13] = x[13] ^ x[1]
vmov s25, r6    @ spill x[25] from r6
vmov r6, s29    @ load x[29] from s29
eor r8, r8, r6    @ x[20] = x[20] ^ x[29]
vmov s29, r6    @ spill x[29] from r6
vmov r6, s3    @ load x[3] from s3
eor r9, r9, r6    @ x[1] = x[1] ^ x[3]
eor r0, r0, r3    @ x[22] = x[22] ^ x[21]
vmov s21, r3    @ spill x[21] from r3
vmov r3, s8    @ load x[8] from s8
eor r3, r3, r5    @ x[8] = x[8] ^ x[9]
vmov s20, r8    @ spill x[20] from r8
vmov r8, s18    @ load x[18] from s18
eor r8, r8, r6    @ x[18] = x[18] ^ x[3]
eor r10, r10, r1    @ x[12] = x[12] ^ x[13]
vmov s12, r10    @ spill x[12] from r10
vmov r10, s17    @ load x[17] from s17
eor r1, r1, r10    @ x[13] = x[13] ^ x[17]
vmov s13, r1    @ spill x[13] from r1
vmov r1, s26    @ load x[26] from s26
eor r1, r1, r9    @ x[26] = x[26] ^ x[1]
eor r9, r9, r5    @ x[1] = x[1] ^ x[9]
vmov s9, r5    @ spill x[9] from r5
vmov r5, s28    @ load x[28] from s28
eor r7, r7, r5    @ x[14] = x[14] ^ x[28]
vmov s14, r7    @ spill x[14] from r7
vmov r7, s23    @ load x[23] from s23
vmov s18, r8    @ spill x[18] from r8
vmov r8, s16    @ load x[16] from s16
eor r7, r7, r8    @ x[23] = x[23] ^ x[16]
vmov s17, r10    @ spill x[17] from r10
vmov r10, s30    @ load x[30] from s30
eor r11, r11, r10    @ x[2] = x[2] ^ x[30]
vmov s26, r1    @ spill x[26] from r1
vmov r1, s10    @ load x[10] from s10
eor r3, r3, r1    @ x[8] = x[8] ^ x[10]
eor r4, r4, r2    @ x[0] = x[0] ^ x[15]
eor r6, r6, r7    @ x[3] = x[3] ^ x[23]
vmov s0, r4    @ spill x[0] from r4
vmov r4, s27    @ load x[27] from s27
vmov s15, r2    @ spill x[15] from r2
vmov r2, s19    @ load x[19] from s19
eor r4, r4, r2    @ x[27] = x[27] ^ x[19]
vmov s8, r3    @ spill x[8] from r3
vmov r3, s31    @ load x[31] from s31
eor r0, r0, r3    @ x[22] = x[22] ^ x[31]
eor r7, r7, r1    @ x[23] = x[23] ^ x[10]
vmov s28, r5    @ spill x[28] from r5
vmov r5, s25    @ load x[25] from s25
eor r9, r9, r5    @ x[1] = x[1] ^ x[25]
vmov s23, r7    @ spill x[23] from r7
vmov r7, s20    @ load x[20] from s20
eor r10, r10, r7    @ x[30] = x[30] ^ x[20]
eor r9, r9, r11    @ x[1] = x[1] ^ x[2]
vmov s2, r11    @ spill x[2] from r11
vmov r11, s29    @ load x[29] from s29
eor r11, r11, r6    @ x[29] = x[29] ^ x[3]
vmov s1, r9    @ spill x[1] from r9
vmov r9, s11    @ load x[11] from s11
eor r9, r9, r0    @ x[11] = x[11] ^ x[22]
vstr.32  s4, [r14, #0]    @ y[0] = x[4] from s4
str  r4, [r14, #4]    @ y[1] = x[27] from r4
vstr.32  s18, [r14, #8]    @ y[2] = x[18] from s18
str  r6, [r14, #12]    @ y[3] = x[3] from r6
vstr.32  s5, [r14, #16]    @ y[4] = x[5] from s5
vstr.32  s24, [r14, #20]    @ y[5] = x[24] from s24
vstr.32  s15, [r14, #24]    @ y[6] = x[15] from s15
vstr.32  s6, [r14, #28]    @ y[7] = x[6] from s6
str  r2, [r14, #32]    @ y[8] = x[19] from r2
vstr.32  s23, [r14, #36]    @ y[9] = x[23] from s23
str  r7, [r14, #40]    @ y[10] = x[20] from r7
vstr.32  s1, [r14, #44]    @ y[11] = x[1] from s1
str  r0, [r14, #48]    @ y[12] = x[22] from r0
str  r5, [r14, #52]    @ y[13] = x[25] from r5
vstr.32  s12, [r14, #56]    @ y[14] = x[12] from s12
vstr.32  s14, [r14, #60]    @ y[15] = x[14] from s14
str  r8, [r14, #64]    @ y[16] = x[16] from r8
vstr.32  s0, [r14, #68]    @ y[17] = x[0] from s0
str  r9, [r14, #72]    @ y[18] = x[11] from r9
vstr.32  s2, [r14, #76]    @ y[19] = x[2] from s2
str  r3, [r14, #80]    @ y[20] = x[31] from r3
vstr.32  s21, [r14, #84]    @ y[21] = x[21] from s21
vstr.32  s26, [r14, #88]    @ y[22] = x[26] from s26
str  r11, [r14, #92]    @ y[23] = x[29] from r11
vstr.32  s28, [r14, #96]    @ y[24] = x[28] from s28
str  r1, [r14, #100]    @ y[25] = x[10] from r1
vstr.32  s8, [r14, #104]    @ y[26] = x[8] from s8
vstr.32  s13, [r14, #108]    @ y[27] = x[13] from s13
vstr.32  s17, [r14, #112]    @ y[28] = x[17] from s17
vstr.32  s7, [r14, #116]    @ y[29] = x[7] from s7
vstr.32  s9, [r14, #120]    @ y[30] = x[9] from s9
str  r10, [r14, #124]    @ y[31] = x[30] from r10
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v24, .-gft_mul_v24
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v25
.type gft_mul_v25, %function
.align 2
gft_mul_v25:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s11    @ load x[11] from s11
vmov r1, s13    @ load x[13] from s13
eor r0, r0, r1    @ x[11] = x[11] ^ x[13]
vmov r2, s30    @ load x[30] from s30
eor r2, r2, r0    @ x[30] = x[30] ^ x[11]
vmov r3, s22    @ load x[22] from s22
vmov r4, s3    @ load x[3] from s3
eor r3, r3, r4    @ x[22] = x[22] ^ x[3]
vmov r5, s6    @ load x[6] from s6
eor r5, r5, r1    @ x[6] = x[6] ^ x[13]
vmov r6, s5    @ load x[5] from s5
vmov r7, s25    @ load x[25] from s25
eor r6, r6, r7    @ x[5] = x[5] ^ x[25]
vmov r8, s4    @ load x[4] from s4
eor r8, r8, r2    @ x[4] = x[4] ^ x[30]
eor r7, r7, r5    @ x[25] = x[25] ^ x[6]
vmov r9, s26    @ load x[26] from s26
vmov r10, s29    @ load x[29] from s29
eor r9, r9, r10    @ x[26] = x[26] ^ x[29]
vmov r11, s10    @ load x[10] from s10
vmov s30, r2    @ spill x[30] from r2
vmov r2, s31    @ load x[31] from s31
eor r11, r11, r2    @ x[10] = x[10] ^ x[31]
eor r10, r10, r3    @ x[29] = x[29] ^ x[22]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s19    @ load x[19] from s19
vmov s10, r11    @ spill x[10] from r11
vmov r11, s14    @ load x[14] from s14
eor r2, r2, r11    @ x[19] = x[19] ^ x[14]
vmov s26, r9    @ spill x[26] from r9
vmov r9, s24    @ load x[24] from s24
eor r9, r9, r2    @ x[24] = x[24] ^ x[19]
vmov s24, r9    @ spill x[24] from r9
vmov r9, s7    @ load x[7] from s7
eor r11, r11, r9    @ x[14] = x[14] ^ x[7]
vmov s14, r11    @ spill x[14] from r11
vmov r11, s9    @ load x[9] from s9
eor r11, r11, r5    @ x[9] = x[9] ^ x[6]
eor r0, r0, r7    @ x[11] = x[11] ^ x[25]
eor r7, r7, r8    @ x[25] = x[25] ^ x[4]
vmov s9, r11    @ spill x[9] from r11
vmov r11, s2    @ load x[2] from s2
eor r11, r11, r8    @ x[2] = x[2] ^ x[4]
vmov s11, r0    @ spill x[11] from r0
vmov r0, s17    @ load x[17] from s17
eor r3, r3, r0    @ x[22] = x[22] ^ x[17]
eor r8, r8, r6    @ x[4] = x[4] ^ x[5]
vmov s6, r5    @ spill x[6] from r5
vmov r5, s15    @ load x[15] from s15
eor r3, r3, r5    @ x[22] = x[22] ^ x[15]
vmov s17, r0    @ spill x[17] from r0
vmov r0, s0    @ load x[0] from s0
eor r4, r4, r0    @ x[3] = x[3] ^ x[0]
vmov s3, r4    @ spill x[3] from r4
vmov r4, s16    @ load x[16] from s16
eor r5, r5, r4    @ x[15] = x[15] ^ x[16]
eor r4, r4, r10    @ x[16] = x[16] ^ x[29]
eor r10, r10, r7    @ x[29] = x[29] ^ x[25]
eor r0, r0, r3    @ x[0] = x[0] ^ x[22]
eor r3, r3, r7    @ x[22] = x[22] ^ x[25]
eor r5, r5, r6    @ x[15] = x[15] ^ x[5]
eor r7, r7, r6    @ x[25] = x[25] ^ x[5]
vmov s0, r0    @ spill x[0] from r0
vmov r0, s20    @ load x[20] from s20
eor r6, r6, r0    @ x[5] = x[5] ^ x[20]
eor r10, r10, r11    @ x[29] = x[29] ^ x[2]
eor r0, r0, r1    @ x[20] = x[20] ^ x[13]
vmov s25, r7    @ spill x[25] from r7
vmov r7, s23    @ load x[23] from s23
eor r11, r11, r7    @ x[2] = x[2] ^ x[23]
vmov s2, r11    @ spill x[2] from r11
vmov r11, s12    @ load x[12] from s12
vmov s15, r5    @ spill x[15] from r5
vmov r5, s26    @ load x[26] from s26
eor r11, r11, r5    @ x[12] = x[12] ^ x[26]
vmov s12, r11    @ spill x[12] from r11
vmov r11, s1    @ load x[1] from s1
eor r10, r10, r11    @ x[29] = x[29] ^ x[1]
eor r4, r4, r9    @ x[16] = x[16] ^ x[7]
eor r2, r2, r11    @ x[19] = x[19] ^ x[1]
eor r1, r1, r7    @ x[13] = x[13] ^ x[23]
vmov s13, r1    @ spill x[13] from r1
vmov r1, s27    @ load x[27] from s27
eor r8, r8, r1    @ x[4] = x[4] ^ x[27]
vmov s1, r11    @ spill x[1] from r11
vmov r11, s17    @ load x[17] from s17
eor r5, r5, r11    @ x[26] = x[26] ^ x[17]
vmov s26, r5    @ spill x[26] from r5
vmov r5, s6    @ load x[6] from s6
vmov s27, r1    @ spill x[27] from r1
vmov r1, s28    @ load x[28] from s28
eor r5, r5, r1    @ x[6] = x[6] ^ x[28]
vmov s6, r5    @ spill x[6] from r5
vmov r5, s8    @ load x[8] from s8
eor r5, r5, r10    @ x[8] = x[8] ^ x[29]
vmov s28, r1    @ spill x[28] from r1
vmov r1, s11    @ load x[11] from s11
eor r1, r1, r5    @ x[11] = x[11] ^ x[8]
eor r0, r0, r11    @ x[20] = x[20] ^ x[17]
vmov s29, r10    @ spill x[29] from r10
vmov r10, s30    @ load x[30] from s30
eor r11, r11, r10    @ x[17] = x[17] ^ x[30]
eor r9, r9, r3    @ x[7] = x[7] ^ x[22]
eor r3, r3, r7    @ x[22] = x[22] ^ x[23]
vmov s22, r3    @ spill x[22] from r3
vmov r3, s10    @ load x[10] from s10
eor r4, r4, r3    @ x[16] = x[16] ^ x[10]
eor r5, r5, r10    @ x[8] = x[8] ^ x[30]
eor r7, r7, r6    @ x[23] = x[23] ^ x[5]
eor r11, r11, r6    @ x[17] = x[17] ^ x[5]
eor r10, r10, r0    @ x[30] = x[30] ^ x[20]
vmov s30, r10    @ spill x[30] from r10
vmov r10, s24    @ load x[24] from s24
eor r3, r3, r10    @ x[10] = x[10] ^ x[24]
eor r2, r2, r8    @ x[19] = x[19] ^ x[4]
eor r10, r10, r1    @ x[24] = x[24] ^ x[11]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s21    @ load x[21] from s21
vmov s8, r5    @ spill x[8] from r5
vmov r5, s29    @ load x[29] from s29
eor r2, r2, r5    @ x[21] = x[21] ^ x[29]
vmov s20, r0    @ spill x[20] from r0
vmov r0, s14    @ load x[14] from s14
eor r6, r6, r0    @ x[5] = x[5] ^ x[14]
vmov s21, r2    @ spill x[21] from r2
vmov r2, s9    @ load x[9] from s9
eor r6, r6, r2    @ x[5] = x[5] ^ x[9]
eor r2, r2, r8    @ x[9] = x[9] ^ x[4]
eor r1, r1, r0    @ x[11] = x[11] ^ x[14]
vmov s9, r2    @ spill x[9] from r2
vmov r2, s28    @ load x[28] from s28
eor r7, r7, r2    @ x[23] = x[23] ^ x[28]
vmov s5, r6    @ spill x[5] from r6
vmov r6, s15    @ load x[15] from s15
eor r8, r8, r6    @ x[4] = x[4] ^ x[15]
eor r2, r2, r6    @ x[28] = x[28] ^ x[15]
vmov s10, r3    @ spill x[10] from r3
vmov r3, s25    @ load x[25] from s25
eor r0, r0, r3    @ x[14] = x[14] ^ x[25]
vmov s14, r0    @ spill x[14] from r0
vmov r0, s2    @ load x[2] from s2
vmov s28, r2    @ spill x[28] from r2
vmov r2, s31    @ load x[31] from s31
eor r0, r0, r2    @ x[2] = x[2] ^ x[31]
eor r3, r3, r5    @ x[25] = x[25] ^ x[29]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s6    @ load x[6] from s6
eor r3, r3, r10    @ x[6] = x[6] ^ x[24]
vmov s4, r8    @ spill x[4] from r8
vmov r8, s27    @ load x[27] from s27
eor r10, r10, r8    @ x[24] = x[24] ^ x[27]
eor r6, r6, r10    @ x[15] = x[15] ^ x[24]
eor r10, r10, r2    @ x[24] = x[24] ^ x[31]
eor r9, r9, r7    @ x[7] = x[7] ^ x[23]
vmov s7, r9    @ spill x[7] from r9
vmov r9, s1    @ load x[1] from s1
eor r8, r8, r9    @ x[27] = x[27] ^ x[1]
vmov s24, r10    @ spill x[24] from r10
vmov r10, s18    @ load x[18] from s18
eor r1, r1, r10    @ x[11] = x[11] ^ x[18]
eor r11, r11, r5    @ x[17] = x[17] ^ x[29]
eor r3, r3, r10    @ x[6] = x[6] ^ x[18]
vmov s18, r10    @ spill x[18] from r10
vmov r10, s12    @ load x[12] from s12
eor r6, r6, r10    @ x[15] = x[15] ^ x[12]
eor r7, r7, r3    @ x[23] = x[23] ^ x[6]
eor r8, r8, r4    @ x[27] = x[27] ^ x[16]
eor r4, r4, r0    @ x[16] = x[16] ^ x[2]
vmov s16, r4    @ spill x[16] from r4
vmov r4, s4    @ load x[4] from s4
eor r4, r4, r3    @ x[4] = x[4] ^ x[6]
eor r5, r5, r10    @ x[29] = x[29] ^ x[12]
eor r10, r10, r1    @ x[12] = x[12] ^ x[11]
eor r1, r1, r9    @ x[11] = x[11] ^ x[1]
vmov s12, r10    @ spill x[12] from r10
vmov r10, s28    @ load x[28] from s28
eor r2, r2, r10    @ x[31] = x[31] ^ x[28]
vmov s15, r6    @ spill x[15] from r6
vmov r6, s0    @ load x[0] from s0
eor r9, r9, r6    @ x[1] = x[1] ^ x[0]
eor r6, r6, r10    @ x[0] = x[0] ^ x[28]
eor r1, r1, r10    @ x[11] = x[11] ^ x[28]
vmov s0, r6    @ spill x[0] from r6
vmov r6, s10    @ load x[10] from s10
eor r10, r10, r6    @ x[28] = x[28] ^ x[10]
vmov s1, r9    @ spill x[1] from r9
vmov r9, s21    @ load x[21] from s21
eor r8, r8, r9    @ x[27] = x[27] ^ x[21]
vmov s27, r8    @ spill x[27] from r8
vmov r8, s24    @ load x[24] from s24
vmov s6, r3    @ spill x[6] from r3
vmov r3, s25    @ load x[25] from s25
eor r8, r8, r3    @ x[24] = x[24] ^ x[25]
vmov s24, r8    @ spill x[24] from r8
vmov r8, s26    @ load x[26] from s26
eor r3, r3, r8    @ x[25] = x[25] ^ x[26]
eor r7, r7, r9    @ x[23] = x[23] ^ x[21]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s20    @ load x[20] from s20
eor r3, r3, r9    @ x[20] = x[20] ^ x[21]
eor r8, r8, r7    @ x[26] = x[26] ^ x[23]
eor r7, r7, r11    @ x[23] = x[23] ^ x[17]
vmov s20, r3    @ spill x[20] from r3
vmov r3, s5    @ load x[5] from s5
eor r1, r1, r3    @ x[11] = x[11] ^ x[5]
vmov s23, r7    @ spill x[23] from r7
vmov r7, s22    @ load x[22] from s22
eor r7, r7, r2    @ x[22] = x[22] ^ x[31]
eor r2, r2, r5    @ x[31] = x[31] ^ x[29]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s18    @ load x[18] from s18
eor r5, r5, r7    @ x[29] = x[29] ^ x[18]
vmov s5, r3    @ spill x[5] from r3
vmov r3, s13    @ load x[13] from s13
vmov s26, r8    @ spill x[26] from r8
vmov r8, s3    @ load x[3] from s3
eor r3, r3, r8    @ x[13] = x[13] ^ x[3]
eor r7, r7, r0    @ x[18] = x[18] ^ x[2]
vmov s3, r8    @ spill x[3] from r8
vmov r8, s14    @ load x[14] from s14
eor r11, r11, r8    @ x[17] = x[17] ^ x[14]
eor r0, r0, r10    @ x[2] = x[2] ^ x[28]
vmov s17, r11    @ spill x[17] from r11
vmov r11, s9    @ load x[9] from s9
eor r7, r7, r11    @ x[18] = x[18] ^ x[9]
vmov s9, r11    @ spill x[9] from r11
vmov r11, s8    @ load x[8] from s8
eor r1, r1, r11    @ x[11] = x[11] ^ x[8]
vmov s8, r11    @ spill x[8] from r11
vmov r11, s30    @ load x[30] from s30
vmov s28, r10    @ spill x[28] from r10
vmov r10, s19    @ load x[19] from s19
eor r11, r11, r10    @ x[30] = x[30] ^ x[19]
eor r7, r7, r1    @ x[18] = x[18] ^ x[11]
eor r1, r1, r9    @ x[11] = x[11] ^ x[21]
eor r9, r9, r0    @ x[21] = x[21] ^ x[2]
vmov s11, r1    @ spill x[11] from r1
vmov r1, s7    @ load x[7] from s7
eor r1, r1, r5    @ x[7] = x[7] ^ x[29]
eor r6, r6, r2    @ x[10] = x[10] ^ x[31]
eor r10, r10, r4    @ x[19] = x[19] ^ x[4]
eor r4, r4, r1    @ x[4] = x[4] ^ x[7]
eor r7, r7, r4    @ x[18] = x[18] ^ x[4]
vmov s7, r1    @ spill x[7] from r1
vmov r1, s6    @ load x[6] from s6
eor r6, r6, r1    @ x[10] = x[10] ^ x[6]
eor r6, r6, r10    @ x[10] = x[10] ^ x[19]
eor r8, r8, r9    @ x[14] = x[14] ^ x[21]
vmov s19, r10    @ spill x[19] from r10
vmov r10, s15    @ load x[15] from s15
eor r9, r9, r10    @ x[21] = x[21] ^ x[15]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s26    @ load x[26] from s26
vmov s21, r9    @ spill x[21] from r9
vmov r9, s24    @ load x[24] from s24
eor r2, r2, r9    @ x[26] = x[26] ^ x[24]
vmov s24, r9    @ spill x[24] from r9
vmov r9, s27    @ load x[27] from s27
eor r8, r8, r9    @ x[14] = x[14] ^ x[27]
vmov s27, r9    @ spill x[27] from r9
vmov r9, s1    @ load x[1] from s1
eor r9, r9, r11    @ x[1] = x[1] ^ x[30]
eor r5, r5, r9    @ x[29] = x[29] ^ x[1]
vmov s14, r8    @ spill x[14] from r8
vmov r8, s28    @ load x[28] from s28
eor r8, r8, r3    @ x[28] = x[28] ^ x[13]
vmov s28, r8    @ spill x[28] from r8
vmov r8, s0    @ load x[0] from s0
eor r8, r8, r6    @ x[0] = x[0] ^ x[10]
eor r8, r8, r0    @ x[0] = x[0] ^ x[2]
vmov s0, r8    @ spill x[0] from r8
vmov r8, s17    @ load x[17] from s17
eor r9, r9, r8    @ x[1] = x[1] ^ x[17]
eor r0, r0, r1    @ x[2] = x[2] ^ x[6]
eor r8, r8, r1    @ x[17] = x[17] ^ x[6]
vmov s1, r9    @ spill x[1] from r9
vmov r9, s11    @ load x[11] from s11
eor r1, r1, r9    @ x[6] = x[6] ^ x[11]
eor r0, r0, r7    @ x[2] = x[2] ^ x[18]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s5    @ load x[5] from s5
eor r9, r9, r1    @ x[11] = x[11] ^ x[5]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s16    @ load x[16] from s16
eor r6, r6, r7    @ x[16] = x[16] ^ x[18]
eor r7, r7, r5    @ x[18] = x[18] ^ x[29]
eor r9, r9, r10    @ x[11] = x[11] ^ x[15]
eor r9, r9, r11    @ x[11] = x[11] ^ x[30]
eor r10, r10, r0    @ x[15] = x[15] ^ x[2]
eor r10, r10, r6    @ x[15] = x[15] ^ x[16]
eor r5, r5, r3    @ x[29] = x[29] ^ x[13]
eor r1, r1, r8    @ x[5] = x[5] ^ x[17]
vmov s15, r10    @ spill x[15] from r10
vmov r10, s8    @ load x[8] from s8
eor r3, r3, r10    @ x[13] = x[13] ^ x[8]
eor r8, r8, r4    @ x[17] = x[17] ^ x[4]
eor r4, r4, r2    @ x[4] = x[4] ^ x[26]
vmov s13, r3    @ spill x[13] from r3
vmov r3, s21    @ load x[21] from s21
eor r2, r2, r3    @ x[26] = x[26] ^ x[21]
vmov s5, r1    @ spill x[5] from r1
vmov r1, s10    @ load x[10] from s10
eor r3, r3, r1    @ x[21] = x[21] ^ x[10]
eor r1, r1, r7    @ x[10] = x[10] ^ x[18]
vmov s21, r3    @ spill x[21] from r3
vmov r3, s3    @ load x[3] from s3
eor r7, r7, r3    @ x[18] = x[18] ^ x[3]
vmov s11, r9    @ spill x[11] from r9
vmov r9, s24    @ load x[24] from s24
eor r3, r3, r9    @ x[3] = x[3] ^ x[24]
eor r9, r9, r6    @ x[24] = x[24] ^ x[16]
eor r9, r9, r10    @ x[24] = x[24] ^ x[8]
eor r6, r6, r7    @ x[16] = x[16] ^ x[18]
vmov s16, r6    @ spill x[16] from r6
vmov r6, s14    @ load x[14] from s14
eor r7, r7, r6    @ x[18] = x[18] ^ x[14]
vmov s24, r9    @ spill x[24] from r9
vmov r9, s31    @ load x[31] from s31
eor r2, r2, r9    @ x[26] = x[26] ^ x[31]
eor r8, r8, r0    @ x[17] = x[17] ^ x[2]
vmov s17, r8    @ spill x[17] from r8
vmov r8, s25    @ load x[25] from s25
eor r10, r10, r8    @ x[8] = x[8] ^ x[25]
eor r1, r1, r11    @ x[10] = x[10] ^ x[30]
vmov s31, r9    @ spill x[31] from r9
vmov r9, s23    @ load x[23] from s23
eor r0, r0, r9    @ x[2] = x[2] ^ x[23]
eor r6, r6, r1    @ x[14] = x[14] ^ x[10]
eor r1, r1, r8    @ x[10] = x[10] ^ x[25]
vmov s10, r1    @ spill x[10] from r1
vmov r1, s27    @ load x[27] from s27
eor r4, r4, r1    @ x[4] = x[4] ^ x[27]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s22    @ load x[22] from s22
eor r10, r10, r6    @ x[8] = x[8] ^ x[22]
eor r7, r7, r5    @ x[18] = x[18] ^ x[29]
eor r1, r1, r11    @ x[27] = x[27] ^ x[30]
eor r0, r0, r6    @ x[2] = x[2] ^ x[22]
vmov s29, r5    @ spill x[29] from r5
vmov r5, s1    @ load x[1] from s1
vmov s18, r7    @ spill x[18] from r7
vmov r7, s20    @ load x[20] from s20
eor r5, r5, r7    @ x[1] = x[1] ^ x[20]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s0    @ load x[0] from s0
eor r11, r11, r6    @ x[0] = x[0] ^ x[22]
vmov s27, r1    @ spill x[27] from r1
vmov r1, s9    @ load x[9] from s9
vmov s25, r8    @ spill x[25] from r8
vmov r8, s19    @ load x[19] from s19
eor r1, r1, r8    @ x[9] = x[9] ^ x[19]
vmov s2, r0    @ spill x[2] from r0
vmov r0, s6    @ load x[6] from s6
eor r7, r7, r0    @ x[20] = x[20] ^ x[6]
eor r3, r3, r0    @ x[3] = x[3] ^ x[6]
vmov s3, r3    @ spill x[3] from r3
vmov r3, s11    @ load x[11] from s11
vmov s6, r0    @ spill x[6] from r0
vmov r0, s28    @ load x[28] from s28
eor r3, r3, r0    @ x[11] = x[11] ^ x[28]
vmov s28, r0    @ spill x[28] from r0
vmov r0, s21    @ load x[21] from s21
vmov s11, r3    @ spill x[11] from r3
vmov r3, s5    @ load x[5] from s5
eor r0, r0, r3    @ x[21] = x[21] ^ x[5]
eor r11, r11, r2    @ x[0] = x[0] ^ x[26]
vmov s0, r11    @ spill x[0] from r11
vmov r11, s24    @ load x[24] from s24
vmov s21, r0    @ spill x[21] from r0
vmov r0, s31    @ load x[31] from s31
eor r11, r11, r0    @ x[24] = x[24] ^ x[31]
eor r11, r11, r4    @ x[24] = x[24] ^ x[4]
eor r2, r2, r11    @ x[26] = x[26] ^ x[24]
vmov s26, r2    @ spill x[26] from r2
vmov r2, s17    @ load x[17] from s17
eor r2, r2, r7    @ x[17] = x[17] ^ x[20]
vmov s24, r11    @ spill x[24] from r11
vmov r11, s7    @ load x[7] from s7
vmov s20, r7    @ spill x[20] from r7
vmov r7, s13    @ load x[13] from s13
eor r11, r11, r7    @ x[7] = x[7] ^ x[13]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s15    @ load x[15] from s15
eor r2, r2, r5    @ x[15] = x[15] ^ x[1]
eor r1, r1, r4    @ x[9] = x[9] ^ x[4]
eor r8, r8, r11    @ x[19] = x[19] ^ x[7]
eor r11, r11, r4    @ x[7] = x[7] ^ x[4]
eor r9, r9, r10    @ x[23] = x[23] ^ x[8]
eor r4, r4, r3    @ x[4] = x[4] ^ x[5]
eor r6, r6, r5    @ x[22] = x[22] ^ x[1]
eor r7, r7, r1    @ x[13] = x[13] ^ x[9]
eor r0, r0, r1    @ x[31] = x[31] ^ x[9]
vmov s4, r4    @ spill x[4] from r4
vmov r4, s17    @ load x[17] from s17
eor r8, r8, r4    @ x[19] = x[19] ^ x[17]
vmov s8, r10    @ spill x[8] from r10
vmov r10, s11    @ load x[11] from s11
vmov s23, r9    @ spill x[23] from r9
vmov r9, s2    @ load x[2] from s2
eor r10, r10, r9    @ x[11] = x[11] ^ x[2]
vmov s22, r6    @ spill x[22] from r6
vmov r6, s12    @ load x[12] from s12
eor r6, r6, r1    @ x[12] = x[12] ^ x[9]
vmov s1, r5    @ spill x[1] from r5
vmov r5, s20    @ load x[20] from s20
eor r5, r5, r11    @ x[20] = x[20] ^ x[7]
vmov s9, r1    @ spill x[9] from r1
vmov r1, s25    @ load x[25] from s25
vmov s19, r8    @ spill x[19] from r8
vmov r8, s14    @ load x[14] from s14
eor r1, r1, r8    @ x[25] = x[25] ^ x[14]
vstr.32  s22, [r14, #0]    @ y[0] = x[22] from s22
vstr.32  s16, [r14, #4]    @ y[1] = x[16] from s16
vstr.32  s26, [r14, #8]    @ y[2] = x[26] from s26
vstr.32  s21, [r14, #12]    @ y[3] = x[21] from s21
str  r11, [r14, #16]    @ y[4] = x[7] from r11
vstr.32  s23, [r14, #20]    @ y[5] = x[23] from s23
str  r6, [r14, #24]    @ y[6] = x[12] from r6
str  r1, [r14, #28]    @ y[7] = x[25] from r1
vstr.32  s3, [r14, #32]    @ y[8] = x[3] from s3
vstr.32  s19, [r14, #36]    @ y[9] = x[19] from s19
str  r8, [r14, #40]    @ y[10] = x[14] from r8
str  r0, [r14, #44]    @ y[11] = x[31] from r0
vstr.32  s28, [r14, #48]    @ y[12] = x[28] from s28
vstr.32  s24, [r14, #52]    @ y[13] = x[24] from s24
str  r5, [r14, #56]    @ y[14] = x[20] from r5
str  r9, [r14, #60]    @ y[15] = x[2] from r9
vstr.32  s1, [r14, #64]    @ y[16] = x[1] from s1
vstr.32  s29, [r14, #68]    @ y[17] = x[29] from s29
str  r2, [r14, #72]    @ y[18] = x[15] from r2
vstr.32  s30, [r14, #76]    @ y[19] = x[30] from s30
vstr.32  s27, [r14, #80]    @ y[20] = x[27] from s27
str  r10, [r14, #84]    @ y[21] = x[11] from r10
str  r7, [r14, #88]    @ y[22] = x[13] from r7
vstr.32  s10, [r14, #92]    @ y[23] = x[10] from s10
vstr.32  s6, [r14, #96]    @ y[24] = x[6] from s6
vstr.32  s0, [r14, #100]    @ y[25] = x[0] from s0
vstr.32  s8, [r14, #104]    @ y[26] = x[8] from s8
vstr.32  s4, [r14, #108]    @ y[27] = x[4] from s4
vstr.32  s9, [r14, #112]    @ y[28] = x[9] from s9
str  r3, [r14, #116]    @ y[29] = x[5] from r3
str  r4, [r14, #120]    @ y[30] = x[17] from r4
vstr.32  s18, [r14, #124]    @ y[31] = x[18] from s18
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v25, .-gft_mul_v25
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v26
.type gft_mul_v26, %function
.align 2
gft_mul_v26:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s2    @ load x[2] from s2
vmov r1, s6    @ load x[6] from s6
eor r0, r0, r1    @ x[2] = x[2] ^ x[6]
vmov r2, s17    @ load x[17] from s17
vmov r3, s26    @ load x[26] from s26
eor r2, r2, r3    @ x[17] = x[17] ^ x[26]
vmov r4, s23    @ load x[23] from s23
eor r4, r4, r2    @ x[23] = x[23] ^ x[17]
vmov r5, s10    @ load x[10] from s10
eor r0, r0, r5    @ x[2] = x[2] ^ x[10]
vmov r6, s14    @ load x[14] from s14
eor r6, r6, r0    @ x[14] = x[14] ^ x[2]
vmov r7, s21    @ load x[21] from s21
vmov r8, s28    @ load x[28] from s28
eor r7, r7, r8    @ x[21] = x[21] ^ x[28]
vmov r9, s22    @ load x[22] from s22
vmov r10, s1    @ load x[1] from s1
eor r9, r9, r10    @ x[22] = x[22] ^ x[1]
vmov r11, s13    @ load x[13] from s13
vmov s26, r3    @ spill x[26] from r3
vmov r3, s29    @ load x[29] from s29
eor r11, r11, r3    @ x[13] = x[13] ^ x[29]
vmov s21, r7    @ spill x[21] from r7
vmov r7, s4    @ load x[4] from s4
eor r5, r5, r7    @ x[10] = x[10] ^ x[4]
vmov s1, r10    @ spill x[1] from r10
vmov r10, s18    @ load x[18] from s18
eor r5, r5, r10    @ x[10] = x[10] ^ x[18]
vmov s10, r5    @ spill x[10] from r5
vmov r5, s30    @ load x[30] from s30
eor r1, r1, r5    @ x[6] = x[6] ^ x[30]
vmov s28, r8    @ spill x[28] from r8
vmov r8, s9    @ load x[9] from s9
eor r8, r8, r5    @ x[9] = x[9] ^ x[30]
vmov s4, r7    @ spill x[4] from r7
vmov r7, s11    @ load x[11] from s11
eor r7, r7, r11    @ x[11] = x[11] ^ x[13]
vmov s11, r7    @ spill x[11] from r7
vmov r7, s20    @ load x[20] from s20
eor r11, r11, r7    @ x[13] = x[13] ^ x[20]
eor r11, r11, r9    @ x[13] = x[13] ^ x[22]
vmov s9, r8    @ spill x[9] from r8
vmov r8, s15    @ load x[15] from s15
eor r1, r1, r8    @ x[6] = x[6] ^ x[15]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s24    @ load x[24] from s24
vmov s18, r10    @ spill x[18] from r10
vmov r10, s31    @ load x[31] from s31
eor r1, r1, r10    @ x[24] = x[24] ^ x[31]
vmov s31, r10    @ spill x[31] from r10
vmov r10, s7    @ load x[7] from s7
eor r3, r3, r10    @ x[29] = x[29] ^ x[7]
eor r9, r9, r1    @ x[22] = x[22] ^ x[24]
eor r7, r7, r10    @ x[20] = x[20] ^ x[7]
vmov s29, r3    @ spill x[29] from r3
vmov r3, s27    @ load x[27] from s27
eor r10, r10, r3    @ x[7] = x[7] ^ x[27]
eor r5, r5, r4    @ x[30] = x[30] ^ x[23]
eor r7, r7, r6    @ x[20] = x[20] ^ x[14]
eor r8, r8, r0    @ x[15] = x[15] ^ x[2]
eor r11, r11, r4    @ x[13] = x[13] ^ x[23]
eor r2, r2, r7    @ x[17] = x[17] ^ x[20]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s18    @ load x[18] from s18
eor r7, r7, r2    @ x[20] = x[20] ^ x[18]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s9    @ load x[9] from s9
eor r0, r0, r11    @ x[2] = x[2] ^ x[9]
vmov s9, r11    @ spill x[9] from r11
vmov r11, s4    @ load x[4] from s4
eor r1, r1, r11    @ x[24] = x[24] ^ x[4]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s28    @ load x[28] from s28
vmov s2, r0    @ spill x[2] from r0
vmov r0, s8    @ load x[8] from s8
eor r6, r6, r0    @ x[28] = x[28] ^ x[8]
vmov s24, r1    @ spill x[24] from r1
vmov r1, s16    @ load x[16] from s16
vmov s18, r2    @ spill x[18] from r2
vmov r2, s25    @ load x[25] from s25
eor r1, r1, r2    @ x[16] = x[16] ^ x[25]
eor r0, r0, r2    @ x[8] = x[8] ^ x[25]
vmov s8, r0    @ spill x[8] from r0
vmov r0, s0    @ load x[0] from s0
eor r6, r6, r0    @ x[28] = x[28] ^ x[0]
vmov s28, r6    @ spill x[28] from r6
vmov r6, s1    @ load x[1] from s1
eor r1, r1, r6    @ x[16] = x[16] ^ x[1]
eor r2, r2, r3    @ x[25] = x[25] ^ x[27]
vmov s25, r2    @ spill x[25] from r2
vmov r2, s19    @ load x[19] from s19
eor r2, r2, r5    @ x[19] = x[19] ^ x[30]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s12    @ load x[12] from s12
eor r11, r11, r2    @ x[4] = x[4] ^ x[12]
eor r0, r0, r11    @ x[0] = x[0] ^ x[4]
eor r11, r11, r4    @ x[4] = x[4] ^ x[23]
eor r4, r4, r3    @ x[23] = x[23] ^ x[27]
vmov s12, r2    @ spill x[12] from r2
vmov r2, s21    @ load x[21] from s21
eor r6, r6, r2    @ x[1] = x[1] ^ x[21]
eor r7, r7, r9    @ x[20] = x[20] ^ x[22]
eor r3, r3, r0    @ x[27] = x[27] ^ x[0]
vmov s27, r3    @ spill x[27] from r3
vmov r3, s31    @ load x[31] from s31
eor r7, r7, r3    @ x[20] = x[20] ^ x[31]
eor r8, r8, r10    @ x[15] = x[15] ^ x[7]
eor r3, r3, r2    @ x[31] = x[31] ^ x[21]
eor r6, r6, r10    @ x[1] = x[1] ^ x[7]
vmov s20, r7    @ spill x[20] from r7
vmov r7, s18    @ load x[18] from s18
vmov s15, r8    @ spill x[15] from r8
vmov r8, s3    @ load x[3] from s3
eor r7, r7, r8    @ x[18] = x[18] ^ x[3]
vmov s18, r7    @ spill x[18] from r7
vmov r7, s8    @ load x[8] from s8
vmov s31, r3    @ spill x[31] from r3
vmov r3, s24    @ load x[24] from s24
eor r7, r7, r3    @ x[8] = x[8] ^ x[24]
vmov s1, r6    @ spill x[1] from r6
vmov r6, s2    @ load x[2] from s2
eor r6, r6, r11    @ x[2] = x[2] ^ x[4]
eor r0, r0, r1    @ x[0] = x[0] ^ x[16]
eor r0, r0, r7    @ x[0] = x[0] ^ x[8]
vmov s2, r6    @ spill x[2] from r6
vmov r6, s14    @ load x[14] from s14
eor r11, r11, r6    @ x[4] = x[4] ^ x[14]
eor r7, r7, r11    @ x[8] = x[8] ^ x[4]
eor r11, r11, r8    @ x[4] = x[4] ^ x[3]
vmov s8, r7    @ spill x[8] from r7
vmov r7, s13    @ load x[13] from s13
eor r8, r8, r7    @ x[3] = x[3] ^ x[13]
vmov s14, r6    @ spill x[14] from r6
vmov r6, s10    @ load x[10] from s10
eor r10, r10, r6    @ x[7] = x[7] ^ x[10]
eor r7, r7, r10    @ x[13] = x[13] ^ x[7]
eor r10, r10, r9    @ x[7] = x[7] ^ x[22]
eor r9, r9, r4    @ x[22] = x[22] ^ x[23]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s11    @ load x[11] from s11
eor r5, r5, r6    @ x[30] = x[30] ^ x[11]
eor r4, r4, r5    @ x[23] = x[23] ^ x[30]
eor r5, r5, r2    @ x[30] = x[30] ^ x[21]
eor r9, r9, r2    @ x[22] = x[22] ^ x[21]
vmov s22, r9    @ spill x[22] from r9
vmov r9, s1    @ load x[1] from s1
eor r9, r9, r3    @ x[1] = x[1] ^ x[24]
eor r4, r4, r8    @ x[23] = x[23] ^ x[3]
vmov s30, r5    @ spill x[30] from r5
vmov r5, s9    @ load x[9] from s9
eor r2, r2, r5    @ x[21] = x[21] ^ x[9]
vmov s21, r2    @ spill x[21] from r2
vmov r2, s31    @ load x[31] from s31
eor r8, r8, r2    @ x[3] = x[3] ^ x[31]
eor r1, r1, r2    @ x[16] = x[16] ^ x[31]
eor r11, r11, r7    @ x[4] = x[4] ^ x[13]
eor r5, r5, r4    @ x[9] = x[9] ^ x[23]
eor r4, r4, r7    @ x[23] = x[23] ^ x[13]
vmov s13, r7    @ spill x[13] from r7
vmov r7, s5    @ load x[5] from s5
eor r7, r7, r3    @ x[5] = x[5] ^ x[24]
vmov s5, r7    @ spill x[5] from r7
vmov r7, s12    @ load x[12] from s12
eor r6, r6, r7    @ x[11] = x[11] ^ x[12]
vmov s3, r8    @ spill x[3] from r8
vmov r8, s15    @ load x[15] from s15
eor r5, r5, r8    @ x[9] = x[9] ^ x[15]
vmov s9, r5    @ spill x[9] from r5
vmov r5, s28    @ load x[28] from s28
eor r5, r5, r4    @ x[28] = x[28] ^ x[23]
vmov s24, r3    @ spill x[24] from r3
vmov r3, s20    @ load x[20] from s20
eor r9, r9, r3    @ x[1] = x[1] ^ x[20]
vmov s1, r9    @ spill x[1] from r9
vmov r9, s17    @ load x[17] from s17
eor r9, r9, r6    @ x[17] = x[17] ^ x[11]
vmov s11, r6    @ spill x[11] from r6
vmov r6, s26    @ load x[26] from s26
eor r6, r6, r10    @ x[26] = x[26] ^ x[7]
eor r10, r10, r8    @ x[7] = x[7] ^ x[15]
vmov s20, r3    @ spill x[20] from r3
vmov r3, s19    @ load x[19] from s19
eor r0, r0, r3    @ x[0] = x[0] ^ x[19]
eor r9, r9, r4    @ x[17] = x[17] ^ x[23]
vmov s0, r0    @ spill x[0] from r0
vmov r0, s29    @ load x[29] from s29
eor r4, r4, r0    @ x[23] = x[23] ^ x[29]
eor r8, r8, r7    @ x[15] = x[15] ^ x[12]
eor r2, r2, r8    @ x[31] = x[31] ^ x[15]
eor r11, r11, r1    @ x[4] = x[4] ^ x[16]
vmov s4, r11    @ spill x[4] from r11
vmov r11, s24    @ load x[24] from s24
eor r11, r11, r5    @ x[24] = x[24] ^ x[28]
vmov s24, r11    @ spill x[24] from r11
vmov r11, s30    @ load x[30] from s30
eor r2, r2, r11    @ x[31] = x[31] ^ x[30]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s10    @ load x[10] from s10
eor r11, r11, r2    @ x[30] = x[30] ^ x[10]
vmov s23, r4    @ spill x[23] from r4
vmov r4, s9    @ load x[9] from s9
eor r5, r5, r4    @ x[28] = x[28] ^ x[9]
vmov s28, r5    @ spill x[28] from r5
vmov r5, s14    @ load x[14] from s14
eor r6, r6, r5    @ x[26] = x[26] ^ x[14]
vmov s9, r4    @ spill x[9] from r4
vmov r4, s18    @ load x[18] from s18
eor r8, r8, r4    @ x[15] = x[15] ^ x[18]
vmov s26, r6    @ spill x[26] from r6
vmov r6, s8    @ load x[8] from s8
eor r8, r8, r6    @ x[15] = x[15] ^ x[8]
eor r4, r4, r5    @ x[18] = x[18] ^ x[14]
vmov s15, r8    @ spill x[15] from r8
vmov r8, s27    @ load x[27] from s27
vmov s8, r6    @ spill x[8] from r6
vmov r6, s21    @ load x[21] from s21
eor r8, r8, r6    @ x[27] = x[27] ^ x[21]
eor r4, r4, r3    @ x[18] = x[18] ^ x[19]
vmov s27, r8    @ spill x[27] from r8
vmov r8, s3    @ load x[3] from s3
eor r10, r10, r8    @ x[7] = x[7] ^ x[3]
eor r10, r10, r3    @ x[7] = x[7] ^ x[19]
eor r7, r7, r6    @ x[12] = x[12] ^ x[21]
eor r6, r6, r5    @ x[21] = x[21] ^ x[14]
eor r6, r6, r0    @ x[21] = x[21] ^ x[29]
eor r5, r5, r1    @ x[14] = x[14] ^ x[16]
vmov s7, r10    @ spill x[7] from r10
vmov r10, s6    @ load x[6] from s6
eor r2, r2, r10    @ x[10] = x[10] ^ x[6]
vmov s14, r5    @ spill x[14] from r5
vmov r5, s20    @ load x[20] from s20
eor r5, r5, r10    @ x[20] = x[20] ^ x[6]
vmov s21, r6    @ spill x[21] from r6
vmov r6, s25    @ load x[25] from s25
eor r10, r10, r6    @ x[6] = x[6] ^ x[25]
vmov s12, r7    @ spill x[12] from r7
vmov r7, s5    @ load x[5] from s5
eor r6, r6, r7    @ x[25] = x[25] ^ x[5]
eor r7, r7, r1    @ x[5] = x[5] ^ x[16]
eor r1, r1, r3    @ x[16] = x[16] ^ x[19]
vmov s16, r1    @ spill x[16] from r1
vmov r1, s11    @ load x[11] from s11
eor r3, r3, r1    @ x[19] = x[19] ^ x[11]
eor r1, r1, r2    @ x[11] = x[11] ^ x[10]
eor r2, r2, r5    @ x[10] = x[10] ^ x[20]
eor r6, r6, r11    @ x[25] = x[25] ^ x[30]
eor r11, r11, r0    @ x[30] = x[30] ^ x[29]
eor r0, r0, r8    @ x[29] = x[29] ^ x[3]
eor r8, r8, r9    @ x[3] = x[3] ^ x[17]
eor r2, r2, r9    @ x[10] = x[10] ^ x[17]
vmov s25, r6    @ spill x[25] from r6
vmov r6, s22    @ load x[22] from s22
eor r10, r10, r6    @ x[6] = x[6] ^ x[22]
vmov s6, r10    @ spill x[6] from r10
vmov r10, s23    @ load x[23] from s23
eor r9, r9, r10    @ x[17] = x[17] ^ x[23]
vmov s29, r0    @ spill x[29] from r0
vmov r0, s8    @ load x[8] from s8
eor r6, r6, r0    @ x[22] = x[22] ^ x[8]
vmov s22, r6    @ spill x[22] from r6
vmov r6, s26    @ load x[26] from s26
eor r0, r0, r6    @ x[8] = x[8] ^ x[26]
eor r10, r10, r0    @ x[23] = x[23] ^ x[8]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s1    @ load x[1] from s1
eor r0, r0, r11    @ x[8] = x[8] ^ x[1]
vmov s8, r0    @ spill x[8] from r0
vmov r0, s15    @ load x[15] from s15
vmov s10, r2    @ spill x[10] from r2
vmov r2, s9    @ load x[9] from s9
eor r0, r0, r2    @ x[15] = x[15] ^ x[9]
eor r0, r0, r9    @ x[15] = x[15] ^ x[17]
eor r9, r9, r7    @ x[17] = x[17] ^ x[5]
eor r7, r7, r2    @ x[5] = x[5] ^ x[9]
eor r8, r8, r3    @ x[3] = x[3] ^ x[19]
vmov s15, r0    @ spill x[15] from r0
vmov r0, s4    @ load x[4] from s4
eor r4, r4, r0    @ x[18] = x[18] ^ x[4]
vmov s18, r4    @ spill x[18] from r4
vmov r4, s12    @ load x[12] from s12
vmov s17, r9    @ spill x[17] from r9
vmov r9, s28    @ load x[28] from s28
eor r4, r4, r9    @ x[12] = x[12] ^ x[28]
eor r4, r4, r8    @ x[12] = x[12] ^ x[3]
vmov s3, r8    @ spill x[3] from r8
vmov r8, s13    @ load x[13] from s13
eor r8, r8, r6    @ x[13] = x[13] ^ x[26]
vmov s13, r8    @ spill x[13] from r8
vmov r8, s27    @ load x[27] from s27
eor r1, r1, r8    @ x[11] = x[11] ^ x[27]
vmov s12, r4    @ spill x[12] from r4
vmov r4, s2    @ load x[2] from s2
eor r2, r2, r4    @ x[9] = x[9] ^ x[2]
eor r2, r2, r5    @ x[9] = x[9] ^ x[20]
eor r5, r5, r11    @ x[20] = x[20] ^ x[1]
eor r11, r11, r7    @ x[1] = x[1] ^ x[5]
vmov s9, r2    @ spill x[9] from r2
vmov r2, s31    @ load x[31] from s31
eor r0, r0, r2    @ x[4] = x[4] ^ x[31]
vmov s4, r0    @ spill x[4] from r0
vmov r0, s24    @ load x[24] from s24
eor r10, r10, r0    @ x[23] = x[23] ^ x[24]
eor r10, r10, r7    @ x[23] = x[23] ^ x[5]
eor r7, r7, r3    @ x[5] = x[5] ^ x[19]
eor r3, r3, r4    @ x[19] = x[19] ^ x[2]
vmov s23, r10    @ spill x[23] from r10
vmov r10, s10    @ load x[10] from s10
eor r10, r10, r6    @ x[10] = x[10] ^ x[26]
eor r6, r6, r5    @ x[26] = x[26] ^ x[20]
vmov s10, r10    @ spill x[10] from r10
vmov r10, s21    @ load x[21] from s21
eor r6, r6, r10    @ x[26] = x[26] ^ x[21]
vmov s26, r6    @ spill x[26] from r6
vmov r6, s17    @ load x[17] from s17
vmov s5, r7    @ spill x[5] from r7
vmov r7, s14    @ load x[14] from s14
eor r6, r6, r7    @ x[17] = x[17] ^ x[14]
eor r9, r9, r11    @ x[28] = x[28] ^ x[1]
vmov s19, r3    @ spill x[19] from r3
vmov r3, s18    @ load x[18] from s18
eor r10, r10, r3    @ x[21] = x[21] ^ x[18]
eor r9, r9, r10    @ x[28] = x[28] ^ x[21]
eor r11, r11, r7    @ x[1] = x[1] ^ x[14]
vmov s1, r11    @ spill x[1] from r11
vmov r11, s30    @ load x[30] from s30
eor r7, r7, r11    @ x[14] = x[14] ^ x[30]
eor r11, r11, r3    @ x[30] = x[30] ^ x[18]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s29    @ load x[29] from s29
eor r3, r3, r9    @ x[18] = x[18] ^ x[29]
eor r7, r7, r1    @ x[14] = x[14] ^ x[11]
vmov s18, r3    @ spill x[18] from r3
vmov r3, s16    @ load x[16] from s16
eor r9, r9, r3    @ x[29] = x[29] ^ x[16]
eor r10, r10, r0    @ x[21] = x[21] ^ x[24]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s7    @ load x[7] from s7
eor r9, r9, r11    @ x[29] = x[29] ^ x[7]
eor r3, r3, r1    @ x[16] = x[16] ^ x[11]
vmov s7, r11    @ spill x[7] from r11
vmov r11, s0    @ load x[0] from s0
eor r1, r1, r11    @ x[11] = x[11] ^ x[0]
vmov s11, r1    @ spill x[11] from r1
vmov r1, s8    @ load x[8] from s8
eor r9, r9, r1    @ x[29] = x[29] ^ x[8]
eor r1, r1, r8    @ x[8] = x[8] ^ x[27]
vmov s0, r11    @ spill x[0] from r11
vmov r11, s6    @ load x[6] from s6
eor r7, r7, r11    @ x[14] = x[14] ^ x[6]
vmov s6, r11    @ spill x[6] from r11
vmov r11, s12    @ load x[12] from s12
eor r0, r0, r11    @ x[24] = x[24] ^ x[12]
vmov s24, r0    @ spill x[24] from r0
vmov r0, s25    @ load x[25] from s25
eor r1, r1, r0    @ x[8] = x[8] ^ x[25]
eor r7, r7, r5    @ x[14] = x[14] ^ x[20]
vmov s25, r0    @ spill x[25] from r0
vmov r0, s15    @ load x[15] from s15
eor r5, r5, r0    @ x[20] = x[20] ^ x[15]
vmov s20, r5    @ spill x[20] from r5
vmov r5, s9    @ load x[9] from s9
eor r9, r9, r5    @ x[29] = x[29] ^ x[9]
eor r10, r10, r4    @ x[21] = x[21] ^ x[2]
eor r6, r6, r8    @ x[17] = x[17] ^ x[27]
vmov s15, r0    @ spill x[15] from r0
vmov r0, s3    @ load x[3] from s3
eor r0, r0, r4    @ x[3] = x[3] ^ x[2]
vmov s3, r0    @ spill x[3] from r0
vmov r0, s13    @ load x[13] from s13
eor r3, r3, r0    @ x[16] = x[16] ^ x[13]
vmov s29, r9    @ spill x[29] from r9
vmov r9, s19    @ load x[19] from s19
eor r2, r2, r9    @ x[31] = x[31] ^ x[19]
eor r9, r9, r0    @ x[19] = x[19] ^ x[13]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s5    @ load x[5] from s5
vmov s14, r7    @ spill x[14] from r7
vmov r7, s0    @ load x[0] from s0
eor r2, r2, r7    @ x[5] = x[5] ^ x[0]
vmov s0, r7    @ spill x[0] from r7
vmov r7, s7    @ load x[7] from s7
eor r8, r8, r7    @ x[27] = x[27] ^ x[7]
vmov s27, r8    @ spill x[27] from r8
vmov r8, s11    @ load x[11] from s11
eor r8, r8, r11    @ x[11] = x[11] ^ x[12]
eor r4, r4, r11    @ x[2] = x[2] ^ x[12]
vmov s2, r4    @ spill x[2] from r4
vmov r4, s30    @ load x[30] from s30
eor r5, r5, r4    @ x[9] = x[9] ^ x[30]
vmov s12, r11    @ spill x[12] from r11
vmov r11, s10    @ load x[10] from s10
eor r11, r11, r6    @ x[10] = x[10] ^ x[17]
vmov s9, r5    @ spill x[9] from r5
vmov r5, s23    @ load x[23] from s23
eor r3, r3, r5    @ x[16] = x[16] ^ x[23]
vmov s16, r3    @ spill x[16] from r3
vmov r3, s28    @ load x[28] from s28
eor r2, r2, r3    @ x[5] = x[5] ^ x[28]
eor r0, r0, r10    @ x[13] = x[13] ^ x[21]
vmov s21, r10    @ spill x[21] from r10
vmov r10, s22    @ load x[22] from s22
vmov s5, r2    @ spill x[5] from r2
vmov r2, s26    @ load x[26] from s26
eor r10, r10, r2    @ x[22] = x[22] ^ x[26]
eor r11, r11, r0    @ x[10] = x[10] ^ x[13]
eor r8, r8, r2    @ x[11] = x[11] ^ x[26]
vmov s13, r0    @ spill x[13] from r0
vmov r0, s18    @ load x[18] from s18
eor r0, r0, r9    @ x[18] = x[18] ^ x[19]
vmov s7, r7    @ spill x[7] from r7
vmov r7, s14    @ load x[14] from s14
eor r1, r1, r7    @ x[8] = x[8] ^ x[14]
vmov s11, r8    @ spill x[11] from r8
vmov r8, s29    @ load x[29] from s29
eor r7, r7, r8    @ x[14] = x[14] ^ x[29]
eor r2, r2, r11    @ x[26] = x[26] ^ x[10]
eor r6, r6, r8    @ x[17] = x[17] ^ x[29]
eor r9, r9, r1    @ x[19] = x[19] ^ x[8]
str  r11, [r14, #0]    @ y[0] = x[10] from r11
str  r1, [r14, #4]    @ y[1] = x[8] from r1
vstr.32  s20, [r14, #8]    @ y[2] = x[20] from s20
str  r2, [r14, #12]    @ y[3] = x[26] from r2
vstr.32  s5, [r14, #16]    @ y[4] = x[5] from s5
vstr.32  s24, [r14, #20]    @ y[5] = x[24] from s24
str  r5, [r14, #24]    @ y[6] = x[23] from r5
vstr.32  s2, [r14, #28]    @ y[7] = x[2] from s2
str  r3, [r14, #32]    @ y[8] = x[28] from r3
vstr.32  s13, [r14, #36]    @ y[9] = x[13] from s13
vstr.32  s21, [r14, #40]    @ y[10] = x[21] from s21
vstr.32  s1, [r14, #44]    @ y[11] = x[1] from s1
str  r0, [r14, #48]    @ y[12] = x[18] from r0
str  r6, [r14, #52]    @ y[13] = x[17] from r6
vstr.32  s27, [r14, #56]    @ y[14] = x[27] from s27
str  r10, [r14, #60]    @ y[15] = x[22] from r10
vstr.32  s6, [r14, #64]    @ y[16] = x[6] from s6
vstr.32  s11, [r14, #68]    @ y[17] = x[11] from s11
vstr.32  s16, [r14, #72]    @ y[18] = x[16] from s16
vstr.32  s7, [r14, #76]    @ y[19] = x[7] from s7
vstr.32  s15, [r14, #80]    @ y[20] = x[15] from s15
vstr.32  s0, [r14, #84]    @ y[21] = x[0] from s0
vstr.32  s4, [r14, #88]    @ y[22] = x[4] from s4
vstr.32  s9, [r14, #92]    @ y[23] = x[9] from s9
str  r7, [r14, #96]    @ y[24] = x[14] from r7
str  r8, [r14, #100]    @ y[25] = x[29] from r8
vstr.32  s12, [r14, #104]    @ y[26] = x[12] from s12
vstr.32  s3, [r14, #108]    @ y[27] = x[3] from s3
vstr.32  s31, [r14, #112]    @ y[28] = x[31] from s31
str  r9, [r14, #116]    @ y[29] = x[19] from r9
str  r4, [r14, #120]    @ y[30] = x[30] from r4
vstr.32  s25, [r14, #124]    @ y[31] = x[25] from s25
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v26, .-gft_mul_v26
.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v27
.type gft_mul_v27, %function
.align 2
gft_mul_v27:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s16    @ load x[16] from s16
vmov r1, s4    @ load x[4] from s4
eor r0, r0, r1    @ x[16] = x[16] ^ x[4]
vmov r2, s6    @ load x[6] from s6
eor r2, r2, r0    @ x[6] = x[6] ^ x[16]
vmov r3, s13    @ load x[13] from s13
vmov r4, s23    @ load x[23] from s23
eor r3, r3, r4    @ x[13] = x[13] ^ x[23]
vmov r5, s25    @ load x[25] from s25
eor r5, r5, r2    @ x[25] = x[25] ^ x[6]
vmov r6, s22    @ load x[22] from s22
vmov r7, s12    @ load x[12] from s12
eor r6, r6, r7    @ x[22] = x[22] ^ x[12]
vmov r8, s19    @ load x[19] from s19
vmov r9, s9    @ load x[9] from s9
eor r8, r8, r9    @ x[19] = x[19] ^ x[9]
vmov r10, s5    @ load x[5] from s5
vmov r11, s26    @ load x[26] from s26
eor r10, r10, r11    @ x[5] = x[5] ^ x[26]
vmov s6, r2    @ spill x[6] from r2
vmov r2, s21    @ load x[21] from s21
vmov s16, r0    @ spill x[16] from r0
vmov r0, s3    @ load x[3] from s3
eor r2, r2, r0    @ x[21] = x[21] ^ x[3]
vmov s3, r0    @ spill x[3] from r0
vmov r0, s29    @ load x[29] from s29
eor r0, r0, r10    @ x[29] = x[29] ^ x[5]
vmov s9, r9    @ spill x[9] from r9
vmov r9, s14    @ load x[14] from s14
eor r0, r0, r9    @ x[29] = x[29] ^ x[14]
eor r4, r4, r8    @ x[23] = x[23] ^ x[19]
vmov s23, r4    @ spill x[23] from r4
vmov r4, s15    @ load x[15] from s15
vmov s22, r6    @ spill x[22] from r6
vmov r6, s28    @ load x[28] from s28
eor r4, r4, r6    @ x[15] = x[15] ^ x[28]
eor r6, r6, r3    @ x[28] = x[28] ^ x[13]
eor r0, r0, r5    @ x[29] = x[29] ^ x[25]
vmov s13, r3    @ spill x[13] from r3
vmov r3, s31    @ load x[31] from s31
eor r5, r5, r3    @ x[25] = x[25] ^ x[31]
vmov s15, r4    @ spill x[15] from r4
vmov r4, s18    @ load x[18] from s18
eor r4, r4, r5    @ x[18] = x[18] ^ x[25]
eor r8, r8, r6    @ x[19] = x[19] ^ x[28]
vmov s28, r6    @ spill x[28] from r6
vmov r6, s30    @ load x[30] from s30
vmov s25, r5    @ spill x[25] from r5
vmov r5, s10    @ load x[10] from s10
eor r6, r6, r5    @ x[30] = x[30] ^ x[10]
eor r6, r6, r8    @ x[30] = x[30] ^ x[19]
eor r6, r6, r1    @ x[30] = x[30] ^ x[4]
eor r1, r1, r10    @ x[4] = x[4] ^ x[5]
eor r2, r2, r9    @ x[21] = x[21] ^ x[14]
eor r7, r7, r4    @ x[12] = x[12] ^ x[18]
vmov s4, r1    @ spill x[4] from r1
vmov r1, s0    @ load x[0] from s0
eor r8, r8, r1    @ x[19] = x[19] ^ x[0]
eor r3, r3, r0    @ x[31] = x[31] ^ x[29]
vmov s29, r0    @ spill x[29] from r0
vmov r0, s7    @ load x[7] from s7
eor r0, r0, r11    @ x[7] = x[7] ^ x[26]
vmov s7, r0    @ spill x[7] from r0
vmov r0, s22    @ load x[22] from s22
eor r1, r1, r0    @ x[0] = x[0] ^ x[22]
vmov s19, r8    @ spill x[19] from r8
vmov r8, s24    @ load x[24] from s24
eor r1, r1, r8    @ x[0] = x[0] ^ x[24]
vmov s26, r11    @ spill x[26] from r11
vmov r11, s11    @ load x[11] from s11
vmov s22, r0    @ spill x[22] from r0
vmov r0, s20    @ load x[20] from s20
eor r11, r11, r0    @ x[11] = x[11] ^ x[20]
eor r10, r10, r1    @ x[5] = x[5] ^ x[0]
vmov s5, r10    @ spill x[5] from r10
vmov r10, s9    @ load x[9] from s9
eor r10, r10, r3    @ x[9] = x[9] ^ x[31]
vmov s31, r3    @ spill x[31] from r3
vmov r3, s25    @ load x[25] from s25
vmov s20, r0    @ spill x[20] from r0
vmov r0, s17    @ load x[17] from s17
eor r3, r3, r0    @ x[25] = x[25] ^ x[17]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s28    @ load x[28] from s28
eor r3, r3, r0    @ x[28] = x[28] ^ x[17]
eor r7, r7, r0    @ x[12] = x[12] ^ x[17]
eor r10, r10, r8    @ x[9] = x[9] ^ x[24]
eor r8, r8, r9    @ x[24] = x[24] ^ x[14]
vmov s12, r7    @ spill x[12] from r7
vmov r7, s16    @ load x[16] from s16
eor r7, r7, r11    @ x[16] = x[16] ^ x[11]
eor r5, r5, r4    @ x[10] = x[10] ^ x[18]
eor r4, r4, r3    @ x[18] = x[18] ^ x[28]
eor r11, r11, r6    @ x[11] = x[11] ^ x[30]
vmov s10, r5    @ spill x[10] from r5
vmov r5, s2    @ load x[2] from s2
eor r5, r5, r0    @ x[2] = x[2] ^ x[17]
vmov s16, r7    @ spill x[16] from r7
vmov r7, s1    @ load x[1] from s1
vmov s18, r4    @ spill x[18] from r4
vmov r4, s27    @ load x[27] from s27
eor r7, r7, r4    @ x[1] = x[1] ^ x[27]
eor r1, r1, r3    @ x[0] = x[0] ^ x[28]
eor r1, r1, r2    @ x[0] = x[0] ^ x[21]
vmov s0, r1    @ spill x[0] from r1
vmov r1, s20    @ load x[20] from s20
eor r0, r0, r1    @ x[17] = x[17] ^ x[20]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s8    @ load x[8] from s8
vmov s21, r2    @ spill x[21] from r2
vmov r2, s15    @ load x[15] from s15
eor r4, r4, r2    @ x[8] = x[8] ^ x[15]
vmov s15, r2    @ spill x[15] from r2
vmov r2, s3    @ load x[3] from s3
eor r0, r0, r2    @ x[17] = x[17] ^ x[3]
vmov s3, r2    @ spill x[3] from r2
vmov r2, s22    @ load x[22] from s22
eor r11, r11, r2    @ x[11] = x[11] ^ x[22]
eor r6, r6, r9    @ x[30] = x[30] ^ x[14]
eor r1, r1, r6    @ x[20] = x[20] ^ x[30]
vmov s11, r11    @ spill x[11] from r11
vmov r11, s13    @ load x[13] from s13
eor r6, r6, r11    @ x[30] = x[30] ^ x[13]
eor r5, r5, r8    @ x[2] = x[2] ^ x[24]
eor r7, r7, r11    @ x[1] = x[1] ^ x[13]
vmov s1, r7    @ spill x[1] from r7
vmov r7, s26    @ load x[26] from s26
eor r5, r5, r7    @ x[2] = x[2] ^ x[26]
eor r0, r0, r2    @ x[17] = x[17] ^ x[22]
vmov s17, r0    @ spill x[17] from r0
vmov r0, s23    @ load x[23] from s23
eor r2, r2, r0    @ x[22] = x[22] ^ x[23]
eor r3, r3, r6    @ x[28] = x[28] ^ x[30]
eor r3, r3, r10    @ x[28] = x[28] ^ x[9]
vmov s20, r1    @ spill x[20] from r1
vmov r1, s19    @ load x[19] from s19
eor r11, r11, r1    @ x[13] = x[13] ^ x[19]
eor r10, r10, r7    @ x[9] = x[9] ^ x[26]
eor r6, r6, r2    @ x[30] = x[30] ^ x[22]
eor r2, r2, r7    @ x[22] = x[22] ^ x[26]
eor r7, r7, r1    @ x[26] = x[26] ^ x[19]
vmov s30, r6    @ spill x[30] from r6
vmov r6, s18    @ load x[18] from s18
eor r1, r1, r6    @ x[19] = x[19] ^ x[18]
eor r6, r6, r0    @ x[18] = x[18] ^ x[23]
eor r2, r2, r4    @ x[22] = x[22] ^ x[8]
vmov s22, r2    @ spill x[22] from r2
vmov r2, s25    @ load x[25] from s25
eor r11, r11, r2    @ x[13] = x[13] ^ x[25]
vmov s26, r7    @ spill x[26] from r7
vmov r7, s5    @ load x[5] from s5
eor r2, r2, r7    @ x[25] = x[25] ^ x[5]
vmov s25, r2    @ spill x[25] from r2
vmov r2, s15    @ load x[15] from s15
eor r9, r9, r2    @ x[14] = x[14] ^ x[15]
eor r7, r7, r9    @ x[5] = x[5] ^ x[14]
vmov s24, r8    @ spill x[24] from r8
vmov r8, s29    @ load x[29] from s29
eor r6, r6, r8    @ x[18] = x[18] ^ x[29]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s21    @ load x[21] from s21
eor r7, r7, r6    @ x[5] = x[5] ^ x[21]
eor r9, r9, r4    @ x[14] = x[14] ^ x[8]
vmov s5, r7    @ spill x[5] from r7
vmov r7, s6    @ load x[6] from s6
eor r4, r4, r7    @ x[8] = x[8] ^ x[6]
eor r10, r10, r11    @ x[9] = x[9] ^ x[13]
vmov s9, r10    @ spill x[9] from r10
vmov r10, s7    @ load x[7] from s7
eor r11, r11, r10    @ x[13] = x[13] ^ x[7]
vmov s14, r9    @ spill x[14] from r9
vmov r9, s12    @ load x[12] from s12
eor r3, r3, r9    @ x[28] = x[28] ^ x[12]
vmov s28, r3    @ spill x[28] from r3
vmov r3, s3    @ load x[3] from s3
eor r3, r3, r8    @ x[3] = x[3] ^ x[29]
vmov s3, r3    @ spill x[3] from r3
vmov r3, s4    @ load x[4] from s4
eor r8, r8, r3    @ x[29] = x[29] ^ x[4]
eor r0, r0, r9    @ x[23] = x[23] ^ x[12]
vmov s23, r0    @ spill x[23] from r0
vmov r0, s31    @ load x[31] from s31
eor r0, r0, r5    @ x[31] = x[31] ^ x[2]
eor r5, r5, r7    @ x[2] = x[2] ^ x[6]
eor r9, r9, r3    @ x[12] = x[12] ^ x[4]
vmov s31, r0    @ spill x[31] from r0
vmov r0, s16    @ load x[16] from s16
eor r1, r1, r0    @ x[19] = x[19] ^ x[16]
vmov s19, r1    @ spill x[19] from r1
vmov r1, s24    @ load x[24] from s24
eor r1, r1, r2    @ x[24] = x[24] ^ x[15]
eor r3, r3, r6    @ x[4] = x[4] ^ x[21]
vmov s12, r9    @ spill x[12] from r9
vmov r9, s20    @ load x[20] from s20
vmov s8, r4    @ spill x[8] from r4
vmov r4, s26    @ load x[26] from s26
eor r9, r9, r4    @ x[20] = x[20] ^ x[26]
eor r4, r4, r6    @ x[26] = x[26] ^ x[21]
vmov s20, r9    @ spill x[20] from r9
vmov r9, s30    @ load x[30] from s30
eor r9, r9, r10    @ x[30] = x[30] ^ x[7]
eor r0, r0, r6    @ x[16] = x[16] ^ x[21]
vmov s16, r0    @ spill x[16] from r0
vmov r0, s27    @ load x[27] from s27
eor r8, r8, r0    @ x[29] = x[29] ^ x[27]
vmov s26, r4    @ spill x[26] from r4
vmov r4, s10    @ load x[10] from s10
vmov s30, r9    @ spill x[30] from r9
vmov r9, s1    @ load x[1] from s1
eor r4, r4, r9    @ x[10] = x[10] ^ x[1]
vmov s10, r4    @ spill x[10] from r4
vmov r4, s11    @ load x[11] from s11
eor r7, r7, r4    @ x[6] = x[6] ^ x[11]
eor r6, r6, r0    @ x[21] = x[21] ^ x[27]
eor r6, r6, r10    @ x[21] = x[21] ^ x[7]
eor r3, r3, r5    @ x[4] = x[4] ^ x[2]
vmov s6, r7    @ spill x[6] from r7
vmov r7, s22    @ load x[22] from s22
eor r9, r9, r7    @ x[1] = x[1] ^ x[22]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s18    @ load x[18] from s18
eor r1, r1, r7    @ x[24] = x[24] ^ x[18]
eor r11, r11, r5    @ x[13] = x[13] ^ x[2]
eor r10, r10, r2    @ x[7] = x[7] ^ x[15]
eor r2, r2, r3    @ x[15] = x[15] ^ x[4]
vmov s18, r7    @ spill x[18] from r7
vmov r7, s8    @ load x[8] from s8
eor r8, r8, r7    @ x[29] = x[29] ^ x[8]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s10    @ load x[10] from s10
eor r7, r7, r11    @ x[8] = x[8] ^ x[10]
eor r0, r0, r1    @ x[27] = x[27] ^ x[24]
vmov s1, r9    @ spill x[1] from r9
vmov r9, s17    @ load x[17] from s17
eor r0, r0, r9    @ x[27] = x[27] ^ x[17]
vmov s17, r9    @ spill x[17] from r9
vmov r9, s12    @ load x[12] from s12
eor r3, r3, r9    @ x[4] = x[4] ^ x[12]
vmov s7, r10    @ spill x[7] from r10
vmov r10, s31    @ load x[31] from s31
eor r9, r9, r10    @ x[12] = x[12] ^ x[31]
vmov s31, r10    @ spill x[31] from r10
vmov r10, s25    @ load x[25] from s25
eor r4, r4, r10    @ x[11] = x[11] ^ x[25]
vmov s15, r2    @ spill x[15] from r2
vmov r2, s30    @ load x[30] from s30
eor r0, r0, r2    @ x[27] = x[27] ^ x[30]
eor r9, r9, r10    @ x[12] = x[12] ^ x[25]
vmov s12, r9    @ spill x[12] from r9
vmov r9, s28    @ load x[28] from s28
eor r9, r9, r11    @ x[28] = x[28] ^ x[10]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s14    @ load x[14] from s14
eor r10, r10, r9    @ x[25] = x[25] ^ x[14]
vmov s30, r2    @ spill x[30] from r2
vmov r2, s9    @ load x[9] from s9
eor r9, r9, r2    @ x[14] = x[14] ^ x[9]
eor r0, r0, r2    @ x[27] = x[27] ^ x[9]
eor r9, r9, r7    @ x[14] = x[14] ^ x[8]
eor r5, r5, r11    @ x[2] = x[2] ^ x[10]
vmov s27, r0    @ spill x[27] from r0
vmov r0, s26    @ load x[26] from s26
eor r9, r9, r0    @ x[14] = x[14] ^ x[26]
eor r2, r2, r8    @ x[9] = x[9] ^ x[29]
eor r0, r0, r6    @ x[26] = x[26] ^ x[21]
vmov s9, r2    @ spill x[9] from r2
vmov r2, s16    @ load x[16] from s16
eor r4, r4, r2    @ x[11] = x[11] ^ x[16]
vmov s11, r4    @ spill x[11] from r4
vmov r4, s23    @ load x[23] from s23
eor r0, r0, r4    @ x[26] = x[26] ^ x[23]
vmov s26, r0    @ spill x[26] from r0
vmov r0, s5    @ load x[5] from s5
eor r8, r8, r0    @ x[29] = x[29] ^ x[5]
eor r3, r3, r1    @ x[4] = x[4] ^ x[24]
vmov s4, r3    @ spill x[4] from r3
vmov r3, s15    @ load x[15] from s15
eor r11, r11, r3    @ x[10] = x[10] ^ x[15]
vmov s2, r5    @ spill x[2] from r5
vmov r5, s20    @ load x[20] from s20
eor r10, r10, r5    @ x[25] = x[25] ^ x[20]
vmov s25, r10    @ spill x[25] from r10
vmov r10, s6    @ load x[6] from s6
eor r4, r4, r10    @ x[23] = x[23] ^ x[6]
vmov s8, r7    @ spill x[8] from r7
vmov r7, s7    @ load x[7] from s7
eor r5, r5, r7    @ x[20] = x[20] ^ x[7]
vmov s29, r8    @ spill x[29] from r8
vmov r8, s19    @ load x[19] from s19
eor r5, r5, r8    @ x[20] = x[20] ^ x[19]
vmov s19, r8    @ spill x[19] from r8
vmov r8, s1    @ load x[1] from s1
eor r7, r7, r8    @ x[7] = x[7] ^ x[1]
eor r9, r9, r1    @ x[14] = x[14] ^ x[24]
eor r11, r11, r4    @ x[10] = x[10] ^ x[23]
eor r1, r1, r10    @ x[24] = x[24] ^ x[6]
vmov s20, r5    @ spill x[20] from r5
vmov r5, s31    @ load x[31] from s31
eor r1, r1, r5    @ x[24] = x[24] ^ x[31]
eor r5, r5, r0    @ x[31] = x[31] ^ x[5]
vmov s14, r9    @ spill x[14] from r9
vmov r9, s30    @ load x[30] from s30
eor r4, r4, r9    @ x[23] = x[23] ^ x[30]
vmov s7, r7    @ spill x[7] from r7
vmov r7, s12    @ load x[12] from s12
eor r7, r7, r4    @ x[12] = x[12] ^ x[23]
eor r1, r1, r8    @ x[24] = x[24] ^ x[1]
vmov s24, r1    @ spill x[24] from r1
vmov r1, s17    @ load x[17] from s17
eor r8, r8, r1    @ x[1] = x[1] ^ x[17]
eor r1, r1, r5    @ x[17] = x[17] ^ x[31]
eor r3, r3, r6    @ x[15] = x[15] ^ x[21]
eor r5, r5, r6    @ x[31] = x[31] ^ x[21]
eor r6, r6, r2    @ x[21] = x[21] ^ x[16]
vmov s12, r7    @ spill x[12] from r7
vmov r7, s0    @ load x[0] from s0
eor r2, r2, r7    @ x[16] = x[16] ^ x[0]
eor r7, r7, r10    @ x[0] = x[0] ^ x[6]
vmov s1, r8    @ spill x[1] from r8
vmov r8, s22    @ load x[22] from s22
eor r10, r10, r8    @ x[6] = x[6] ^ x[22]
vmov s31, r5    @ spill x[31] from r5
vmov r5, s28    @ load x[28] from s28
eor r10, r10, r5    @ x[6] = x[6] ^ x[28]
vmov s6, r10    @ spill x[6] from r10
vmov r10, s13    @ load x[13] from s13
eor r1, r1, r10    @ x[17] = x[17] ^ x[13]
vmov s28, r5    @ spill x[28] from r5
vmov r5, s29    @ load x[29] from s29
vmov s13, r10    @ spill x[13] from r10
vmov r10, s3    @ load x[3] from s3
eor r5, r5, r10    @ x[29] = x[29] ^ x[3]
eor r9, r9, r1    @ x[30] = x[30] ^ x[17]
eor r1, r1, r6    @ x[17] = x[17] ^ x[21]
eor r6, r6, r10    @ x[21] = x[21] ^ x[3]
vmov s30, r9    @ spill x[30] from r9
vmov r9, s8    @ load x[8] from s8
eor r0, r0, r9    @ x[5] = x[5] ^ x[8]
vmov s8, r9    @ spill x[8] from r9
vmov r9, s18    @ load x[18] from s18
eor r11, r11, r9    @ x[10] = x[10] ^ x[18]
eor r4, r4, r8    @ x[23] = x[23] ^ x[22]
eor r1, r1, r0    @ x[17] = x[17] ^ x[5]
vmov s5, r0    @ spill x[5] from r0
vmov r0, s7    @ load x[7] from s7
eor r0, r0, r10    @ x[7] = x[7] ^ x[3]
vmov s17, r1    @ spill x[17] from r1
vmov r1, s14    @ load x[14] from s14
eor r1, r1, r0    @ x[14] = x[14] ^ x[7]
eor r0, r0, r7    @ x[7] = x[7] ^ x[0]
vmov s3, r10    @ spill x[3] from r10
vmov r10, s2    @ load x[2] from s2
eor r6, r6, r10    @ x[21] = x[21] ^ x[2]
vmov s18, r9    @ spill x[18] from r9
vmov r9, s27    @ load x[27] from s27
eor r9, r9, r2    @ x[27] = x[27] ^ x[16]
eor r5, r5, r2    @ x[29] = x[29] ^ x[16]
eor r2, r2, r3    @ x[16] = x[16] ^ x[15]
eor r0, r0, r4    @ x[7] = x[7] ^ x[23]
vmov s23, r4    @ spill x[23] from r4
vmov r4, s31    @ load x[31] from s31
eor r4, r4, r1    @ x[31] = x[31] ^ x[14]
eor r2, r2, r4    @ x[16] = x[16] ^ x[31]
vmov s16, r2    @ spill x[16] from r2
vmov r2, s20    @ load x[20] from s20
eor r2, r2, r10    @ x[20] = x[20] ^ x[2]
vmov s15, r3    @ spill x[15] from r3
vmov r3, s19    @ load x[19] from s19
eor r3, r3, r7    @ x[19] = x[19] ^ x[0]
vmov s7, r0    @ spill x[7] from r0
vmov r0, s4    @ load x[4] from s4
eor r8, r8, r0    @ x[22] = x[22] ^ x[4]
vmov s14, r1    @ spill x[14] from r1
vmov r1, s1    @ load x[1] from s1
eor r0, r0, r1    @ x[4] = x[4] ^ x[1]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s9    @ load x[9] from s9
eor r10, r10, r1    @ x[9] = x[9] ^ x[1]
eor r6, r6, r11    @ x[21] = x[21] ^ x[10]
eor r4, r4, r5    @ x[31] = x[31] ^ x[29]
eor r8, r8, r6    @ x[22] = x[22] ^ x[21]
vmov s21, r6    @ spill x[21] from r6
vmov r6, s18    @ load x[18] from s18
eor r6, r6, r5    @ x[18] = x[18] ^ x[29]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s8    @ load x[8] from s8
vmov s22, r8    @ spill x[22] from r8
vmov r8, s12    @ load x[12] from s12
eor r11, r11, r8    @ x[8] = x[8] ^ x[12]
vmov s29, r5    @ spill x[29] from r5
vmov r5, s13    @ load x[13] from s13
eor r5, r5, r9    @ x[13] = x[13] ^ x[27]
eor r3, r3, r11    @ x[19] = x[19] ^ x[8]
eor r7, r7, r2    @ x[0] = x[0] ^ x[20]
eor r1, r1, r5    @ x[1] = x[1] ^ x[13]
str  r2, [r14, #0]    @ y[0] = x[20] from r2
vstr.32  s17, [r14, #4]    @ y[1] = x[17] from s17
str  r11, [r14, #8]    @ y[2] = x[8] from r11
vstr.32  s11, [r14, #12]    @ y[3] = x[11] from s11
str  r7, [r14, #16]    @ y[4] = x[0] from r7
str  r8, [r14, #20]    @ y[5] = x[12] from r8
vstr.32  s24, [r14, #24]    @ y[6] = x[24] from s24
vstr.32  s5, [r14, #28]    @ y[7] = x[5] from s5
vstr.32  s30, [r14, #32]    @ y[8] = x[30] from s30
vstr.32  s23, [r14, #36]    @ y[9] = x[23] from s23
vstr.32  s14, [r14, #40]    @ y[10] = x[14] from s14
vstr.32  s29, [r14, #44]    @ y[11] = x[29] from s29
vstr.32  s6, [r14, #48]    @ y[12] = x[6] from s6
vstr.32  s7, [r14, #52]    @ y[13] = x[7] from s7
str  r9, [r14, #56]    @ y[14] = x[27] from r9
vstr.32  s22, [r14, #60]    @ y[15] = x[22] from s22
vstr.32  s26, [r14, #64]    @ y[16] = x[26] from s26
vstr.32  s28, [r14, #68]    @ y[17] = x[28] from s28
vstr.32  s15, [r14, #72]    @ y[18] = x[15] from s15
vstr.32  s3, [r14, #76]    @ y[19] = x[3] from s3
str  r3, [r14, #80]    @ y[20] = x[19] from r3
str  r4, [r14, #84]    @ y[21] = x[31] from r4
str  r0, [r14, #88]    @ y[22] = x[4] from r0
str  r10, [r14, #92]    @ y[23] = x[9] from r10
str  r5, [r14, #96]    @ y[24] = x[13] from r5
vstr.32  s2, [r14, #100]    @ y[25] = x[2] from s2
vstr.32  s10, [r14, #104]    @ y[26] = x[10] from s10
str  r1, [r14, #108]    @ y[27] = x[1] from r1
vstr.32  s16, [r14, #112]    @ y[28] = x[16] from s16
vstr.32  s25, [r14, #116]    @ y[29] = x[25] from s25
vstr.32  s21, [r14, #120]    @ y[30] = x[21] from s21
str  r6, [r14, #124]    @ y[31] = x[18] from r6
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v27, .-gft_mul_v27

.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_v28
.type gft_mul_v28, %function
.align 2
gft_mul_v28:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldr.32 s0, [r12, #0]
vldr.32 s1, [r12, #4]
vldr.32 s2, [r12, #8]
vldr.32 s3, [r12, #12]
vldr.32 s4, [r12, #16]
vldr.32 s5, [r12, #20]
vldr.32 s6, [r12, #24]
vldr.32 s7, [r12, #28]
vldr.32 s8, [r12, #32]
vldr.32 s9, [r12, #36]
vldr.32 s10, [r12, #40]
vldr.32 s11, [r12, #44]
vldr.32 s12, [r12, #48]
vldr.32 s13, [r12, #52]
vldr.32 s14, [r12, #56]
vldr.32 s15, [r12, #60]
vldr.32 s16, [r12, #64]
vldr.32 s17, [r12, #68]
vldr.32 s18, [r12, #72]
vldr.32 s19, [r12, #76]
vldr.32 s20, [r12, #80]
vldr.32 s21, [r12, #84]
vldr.32 s22, [r12, #88]
vldr.32 s23, [r12, #92]
vldr.32 s24, [r12, #96]
vldr.32 s25, [r12, #100]
vldr.32 s26, [r12, #104]
vldr.32 s27, [r12, #108]
vldr.32 s28, [r12, #112]
vldr.32 s29, [r12, #116]
vldr.32 s30, [r12, #120]
vldr.32 s31, [r12, #124]
vmov r0, s0    @ load x[0] from s0
vmov r1, s21    @ load x[21] from s21
eor r0, r0, r1    @ x[0] = x[0] ^ x[21]
vmov r2, s6    @ load x[6] from s6
vmov r3, s20    @ load x[20] from s20
eor r2, r2, r3    @ x[6] = x[6] ^ x[20]
vmov r4, s24    @ load x[24] from s24
vmov r5, s31    @ load x[31] from s31
eor r4, r4, r5    @ x[24] = x[24] ^ x[31]
vmov r6, s29    @ load x[29] from s29
eor r4, r4, r6    @ x[24] = x[24] ^ x[29]
vmov r7, s27    @ load x[27] from s27
vmov r8, s22    @ load x[22] from s22
eor r7, r7, r8    @ x[27] = x[27] ^ x[22]
vmov r9, s4    @ load x[4] from s4
vmov r10, s2    @ load x[2] from s2
eor r9, r9, r10    @ x[4] = x[4] ^ x[2]
vmov r11, s13    @ load x[13] from s13
vmov s2, r10    @ spill x[2] from r10
vmov r10, s3    @ load x[3] from s3
eor r11, r11, r10    @ x[13] = x[13] ^ x[3]
eor r9, r9, r11    @ x[4] = x[4] ^ x[13]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s30    @ load x[30] from s30
eor r9, r9, r11    @ x[4] = x[4] ^ x[30]
vmov s29, r6    @ spill x[29] from r6
vmov r6, s18    @ load x[18] from s18
vmov s21, r1    @ spill x[21] from r1
vmov r1, s5    @ load x[5] from s5
eor r6, r6, r1    @ x[18] = x[18] ^ x[5]
eor r8, r8, r9    @ x[22] = x[22] ^ x[4]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s1    @ load x[1] from s1
eor r1, r1, r4    @ x[5] = x[5] ^ x[1]
vmov s22, r8    @ spill x[22] from r8
vmov r8, s12    @ load x[12] from s12
eor r8, r8, r2    @ x[12] = x[12] ^ x[6]
vmov s6, r2    @ spill x[6] from r2
vmov r2, s17    @ load x[17] from s17
vmov s30, r11    @ spill x[30] from r11
vmov r11, s10    @ load x[10] from s10
eor r2, r2, r11    @ x[17] = x[17] ^ x[10]
eor r11, r11, r0    @ x[10] = x[10] ^ x[0]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s28    @ load x[28] from s28
eor r2, r2, r11    @ x[17] = x[17] ^ x[28]
eor r8, r8, r0    @ x[12] = x[12] ^ x[0]
vmov s0, r0    @ spill x[0] from r0
vmov r0, s16    @ load x[16] from s16
eor r3, r3, r0    @ x[20] = x[20] ^ x[16]
vmov s20, r3    @ spill x[20] from r3
vmov r3, s9    @ load x[9] from s9
eor r3, r3, r7    @ x[9] = x[9] ^ x[27]
vmov s9, r3    @ spill x[9] from r3
vmov r3, s19    @ load x[19] from s19
eor r1, r1, r3    @ x[5] = x[5] ^ x[19]
eor r0, r0, r6    @ x[16] = x[16] ^ x[18]
eor r3, r3, r7    @ x[19] = x[19] ^ x[27]
vmov s19, r3    @ spill x[19] from r3
vmov r3, s11    @ load x[11] from s11
eor r11, r11, r3    @ x[28] = x[28] ^ x[11]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s8    @ load x[8] from s8
eor r6, r6, r5    @ x[8] = x[8] ^ x[31]
vmov s8, r6    @ spill x[8] from r6
vmov r6, s7    @ load x[7] from s7
eor r4, r4, r6    @ x[1] = x[1] ^ x[7]
eor r11, r11, r4    @ x[28] = x[28] ^ x[1]
eor r5, r5, r9    @ x[31] = x[31] ^ x[4]
eor r9, r9, r10    @ x[4] = x[4] ^ x[3]
eor r10, r10, r8    @ x[3] = x[3] ^ x[12]
vmov s11, r3    @ spill x[11] from r3
vmov r3, s25    @ load x[25] from s25
eor r7, r7, r3    @ x[27] = x[27] ^ x[25]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s30    @ load x[30] from s30
eor r9, r9, r3    @ x[4] = x[4] ^ x[30]
vmov s4, r9    @ spill x[4] from r9
vmov r9, s22    @ load x[22] from s22
eor r0, r0, r9    @ x[16] = x[16] ^ x[22]
vmov s16, r0    @ spill x[16] from r0
vmov r0, s15    @ load x[15] from s15
eor r8, r8, r0    @ x[12] = x[12] ^ x[15]
vmov s12, r8    @ spill x[12] from r8
vmov r8, s10    @ load x[10] from s10
eor r5, r5, r8    @ x[31] = x[31] ^ x[10]
eor r10, r10, r2    @ x[3] = x[3] ^ x[17]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s24    @ load x[24] from s24
eor r1, r1, r2    @ x[5] = x[5] ^ x[24]
vmov s5, r1    @ spill x[5] from r1
vmov r1, s6    @ load x[6] from s6
eor r4, r4, r1    @ x[1] = x[1] ^ x[6]
vmov s1, r4    @ spill x[1] from r4
vmov r4, s0    @ load x[0] from s0
vmov s10, r8    @ spill x[10] from r8
vmov r8, s23    @ load x[23] from s23
eor r4, r4, r8    @ x[0] = x[0] ^ x[23]
eor r9, r9, r4    @ x[22] = x[22] ^ x[0]
vmov s23, r8    @ spill x[23] from r8
vmov r8, s26    @ load x[26] from s26
eor r1, r1, r8    @ x[6] = x[6] ^ x[26]
eor r5, r5, r1    @ x[31] = x[31] ^ x[6]
vmov s26, r8    @ spill x[26] from r8
vmov r8, s21    @ load x[21] from s21
eor r8, r8, r11    @ x[21] = x[21] ^ x[28]
eor r0, r0, r6    @ x[15] = x[15] ^ x[7]
vmov s15, r0    @ spill x[15] from r0
vmov r0, s9    @ load x[9] from s9
eor r1, r1, r0    @ x[6] = x[6] ^ x[9]
eor r7, r7, r2    @ x[27] = x[27] ^ x[24]
vmov s21, r8    @ spill x[21] from r8
vmov r8, s11    @ load x[11] from s11
eor r10, r10, r8    @ x[3] = x[3] ^ x[11]
eor r8, r8, r2    @ x[11] = x[11] ^ x[24]
vmov s3, r10    @ spill x[3] from r10
vmov r10, s14    @ load x[14] from s14
eor r3, r3, r10    @ x[30] = x[30] ^ x[14]
vmov s11, r8    @ spill x[11] from r8
vmov r8, s18    @ load x[18] from s18
eor r8, r8, r10    @ x[18] = x[18] ^ x[14]
vmov s30, r3    @ spill x[30] from r3
vmov r3, s12    @ load x[12] from s12
eor r10, r10, r3    @ x[14] = x[14] ^ x[12]
vmov s9, r0    @ spill x[9] from r0
vmov r0, s4    @ load x[4] from s4
eor r0, r0, r10    @ x[4] = x[4] ^ x[14]
eor r3, r3, r4    @ x[12] = x[12] ^ x[0]
vmov s4, r0    @ spill x[4] from r0
vmov r0, s25    @ load x[25] from s25
eor r0, r0, r6    @ x[25] = x[25] ^ x[7]
eor r4, r4, r5    @ x[0] = x[0] ^ x[31]
vmov s25, r0    @ spill x[25] from r0
vmov r0, s29    @ load x[29] from s29
eor r3, r3, r0    @ x[12] = x[12] ^ x[29]
vmov s14, r10    @ spill x[14] from r10
vmov r10, s10    @ load x[10] from s10
vmov s31, r5    @ spill x[31] from r5
vmov r5, s5    @ load x[5] from s5
eor r10, r10, r5    @ x[10] = x[10] ^ x[5]
vmov s5, r5    @ spill x[5] from r5
vmov r5, s16    @ load x[16] from s16
eor r2, r2, r5    @ x[24] = x[24] ^ x[16]
vmov s16, r5    @ spill x[16] from r5
vmov r5, s1    @ load x[1] from s1
eor r3, r3, r5    @ x[12] = x[12] ^ x[1]
vmov s1, r5    @ spill x[1] from r5
vmov r5, s8    @ load x[8] from s8
eor r3, r3, r5    @ x[12] = x[12] ^ x[8]
vmov s8, r5    @ spill x[8] from r5
vmov r5, s20    @ load x[20] from s20
eor r5, r5, r11    @ x[20] = x[20] ^ x[28]
eor r3, r3, r6    @ x[12] = x[12] ^ x[7]
vmov s12, r3    @ spill x[12] from r3
vmov r3, s19    @ load x[19] from s19
eor r5, r5, r3    @ x[20] = x[20] ^ x[19]
eor r3, r3, r2    @ x[19] = x[19] ^ x[24]
eor r10, r10, r6    @ x[10] = x[10] ^ x[7]
eor r6, r6, r2    @ x[7] = x[7] ^ x[24]
eor r2, r2, r9    @ x[24] = x[24] ^ x[22]
eor r4, r4, r7    @ x[0] = x[0] ^ x[27]
eor r9, r9, r0    @ x[22] = x[22] ^ x[29]
vmov s10, r10    @ spill x[10] from r10
vmov r10, s2    @ load x[2] from s2
eor r7, r7, r10    @ x[27] = x[27] ^ x[2]
eor r10, r10, r1    @ x[2] = x[2] ^ x[6]
eor r1, r1, r8    @ x[6] = x[6] ^ x[18]
eor r8, r8, r0    @ x[18] = x[18] ^ x[29]
vmov s27, r7    @ spill x[27] from r7
vmov r7, s17    @ load x[17] from s17
eor r0, r0, r7    @ x[29] = x[29] ^ x[17]
vmov s0, r4    @ spill x[0] from r4
vmov r4, s8    @ load x[8] from s8
eor r7, r7, r4    @ x[17] = x[17] ^ x[8]
eor r4, r4, r6    @ x[8] = x[8] ^ x[7]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s9    @ load x[9] from s9
eor r6, r6, r10    @ x[7] = x[7] ^ x[9]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s30    @ load x[30] from s30
eor r1, r1, r6    @ x[6] = x[6] ^ x[30]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s5    @ load x[5] from s5
eor r10, r10, r1    @ x[9] = x[9] ^ x[5]
vmov s17, r7    @ spill x[17] from r7
vmov r7, s16    @ load x[16] from s16
eor r6, r6, r7    @ x[30] = x[30] ^ x[16]
vmov s18, r8    @ spill x[18] from r8
vmov r8, s31    @ load x[31] from s31
eor r1, r1, r8    @ x[5] = x[5] ^ x[31]
vmov s5, r1    @ spill x[5] from r1
vmov r1, s11    @ load x[11] from s11
eor r4, r4, r1    @ x[8] = x[8] ^ x[11]
vmov s8, r4    @ spill x[8] from r4
vmov r4, s21    @ load x[21] from s21
eor r7, r7, r4    @ x[16] = x[16] ^ x[21]
vmov s16, r7    @ spill x[16] from r7
vmov r7, s13    @ load x[13] from s13
vmov s9, r10    @ spill x[9] from r10
vmov r10, s26    @ load x[26] from s26
eor r7, r7, r10    @ x[13] = x[13] ^ x[26]
eor r3, r3, r10    @ x[19] = x[19] ^ x[26]
eor r7, r7, r1    @ x[13] = x[13] ^ x[11]
vmov s13, r7    @ spill x[13] from r7
vmov r7, s15    @ load x[15] from s15
eor r0, r0, r7    @ x[29] = x[29] ^ x[15]
eor r4, r4, r2    @ x[21] = x[21] ^ x[24]
vmov s21, r4    @ spill x[21] from r4
vmov r4, s23    @ load x[23] from s23
eor r2, r2, r4    @ x[24] = x[24] ^ x[23]
vmov s24, r2    @ spill x[24] from r2
vmov r2, s14    @ load x[14] from s14
eor r2, r2, r1    @ x[14] = x[14] ^ x[11]
vmov s14, r2    @ spill x[14] from r2
vmov r2, s3    @ load x[3] from s3
eor r4, r4, r2    @ x[23] = x[23] ^ x[3]
eor r2, r2, r6    @ x[3] = x[3] ^ x[30]
eor r6, r6, r8    @ x[30] = x[30] ^ x[31]
eor r1, r1, r4    @ x[11] = x[11] ^ x[23]
eor r4, r4, r10    @ x[23] = x[23] ^ x[26]
eor r6, r6, r9    @ x[30] = x[30] ^ x[22]
eor r5, r5, r3    @ x[20] = x[20] ^ x[19]
vmov s3, r2    @ spill x[3] from r2
vmov r2, s4    @ load x[4] from s4
eor r11, r11, r2    @ x[28] = x[28] ^ x[4]
eor r7, r7, r3    @ x[15] = x[15] ^ x[19]
eor r3, r3, r2    @ x[19] = x[19] ^ x[4]
eor r4, r4, r3    @ x[23] = x[23] ^ x[19]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s18    @ load x[18] from s18
vmov s30, r6    @ spill x[30] from r6
vmov r6, s9    @ load x[9] from s9
eor r7, r7, r6    @ x[18] = x[18] ^ x[9]
vmov s11, r1    @ spill x[11] from r1
vmov r1, s17    @ load x[17] from s17
eor r1, r1, r9    @ x[17] = x[17] ^ x[22]
eor r6, r6, r11    @ x[9] = x[9] ^ x[28]
vmov s28, r11    @ spill x[28] from r11
vmov r11, s1    @ load x[1] from s1
vmov s23, r4    @ spill x[23] from r4
vmov r4, s7    @ load x[7] from s7
eor r11, r11, r4    @ x[1] = x[1] ^ x[7]
vmov s26, r10    @ spill x[26] from r10
vmov r10, s2    @ load x[2] from s2
vmov s17, r1    @ spill x[17] from r1
vmov r1, s12    @ load x[12] from s12
eor r10, r10, r1    @ x[2] = x[2] ^ x[12]
vmov s12, r1    @ spill x[12] from r1
vmov r1, s8    @ load x[8] from s8
vmov s18, r7    @ spill x[18] from r7
vmov r7, s0    @ load x[0] from s0
eor r1, r1, r7    @ x[8] = x[8] ^ x[0]
eor r4, r4, r0    @ x[7] = x[7] ^ x[29]
vmov s7, r4    @ spill x[7] from r4
vmov r4, s24    @ load x[24] from s24
vmov s8, r1    @ spill x[8] from r1
vmov r1, s6    @ load x[6] from s6
eor r4, r4, r1    @ x[24] = x[24] ^ x[6]
eor r7, r7, r5    @ x[0] = x[0] ^ x[20]
eor r3, r3, r11    @ x[19] = x[19] ^ x[1]
eor r11, r11, r9    @ x[1] = x[1] ^ x[22]
eor r9, r9, r2    @ x[22] = x[22] ^ x[4]
eor r7, r7, r10    @ x[0] = x[0] ^ x[2]
eor r3, r3, r6    @ x[19] = x[19] ^ x[9]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s21    @ load x[21] from s21
eor r6, r6, r10    @ x[9] = x[9] ^ x[21]
eor r4, r4, r9    @ x[24] = x[24] ^ x[22]
vmov s19, r3    @ spill x[19] from r3
vmov r3, s27    @ load x[27] from s27
eor r10, r10, r3    @ x[21] = x[21] ^ x[27]
eor r3, r3, r8    @ x[27] = x[27] ^ x[31]
vmov s27, r3    @ spill x[27] from r3
vmov r3, s14    @ load x[14] from s14
eor r8, r8, r3    @ x[31] = x[31] ^ x[14]
eor r7, r7, r3    @ x[0] = x[0] ^ x[14]
vmov s9, r6    @ spill x[9] from r6
vmov r6, s18    @ load x[18] from s18
eor r8, r8, r6    @ x[31] = x[31] ^ x[18]
eor r8, r8, r11    @ x[31] = x[31] ^ x[1]
vmov s0, r7    @ spill x[0] from r7
vmov r7, s13    @ load x[13] from s13
eor r8, r8, r7    @ x[31] = x[31] ^ x[13]
eor r9, r9, r5    @ x[22] = x[22] ^ x[20]
vmov s31, r8    @ spill x[31] from r8
vmov r8, s17    @ load x[17] from s17
eor r10, r10, r8    @ x[21] = x[21] ^ x[17]
eor r9, r9, r6    @ x[22] = x[22] ^ x[18]
vmov s22, r9    @ spill x[22] from r9
vmov r9, s25    @ load x[25] from s25
eor r3, r3, r9    @ x[14] = x[14] ^ x[25]
eor r5, r5, r9    @ x[20] = x[20] ^ x[25]
vmov s20, r5    @ spill x[20] from r5
vmov r5, s26    @ load x[26] from s26
eor r10, r10, r5    @ x[21] = x[21] ^ x[26]
eor r5, r5, r2    @ x[26] = x[26] ^ x[4]
eor r5, r5, r9    @ x[26] = x[26] ^ x[25]
vmov s26, r5    @ spill x[26] from r5
vmov r5, s16    @ load x[16] from s16
eor r2, r2, r5    @ x[4] = x[4] ^ x[16]
vmov s21, r10    @ spill x[21] from r10
vmov r10, s23    @ load x[23] from s23
eor r9, r9, r10    @ x[25] = x[25] ^ x[23]
eor r11, r11, r5    @ x[1] = x[1] ^ x[16]
vmov s1, r11    @ spill x[1] from r11
vmov r11, s11    @ load x[11] from s11
eor r5, r5, r11    @ x[16] = x[16] ^ x[11]
vmov s16, r5    @ spill x[16] from r5
vmov r5, s28    @ load x[28] from s28
eor r11, r11, r5    @ x[11] = x[11] ^ x[28]
eor r10, r10, r6    @ x[23] = x[23] ^ x[18]
eor r6, r6, r5    @ x[18] = x[18] ^ x[28]
eor r2, r2, r0    @ x[4] = x[4] ^ x[29]
eor r2, r2, r1    @ x[4] = x[4] ^ x[6]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s30    @ load x[30] from s30
eor r5, r5, r6    @ x[28] = x[28] ^ x[30]
vmov s4, r2    @ spill x[4] from r2
vmov r2, s5    @ load x[5] from s5
eor r6, r6, r2    @ x[30] = x[30] ^ x[5]
eor r2, r2, r8    @ x[5] = x[5] ^ x[17]
eor r4, r4, r3    @ x[24] = x[24] ^ x[14]
eor r3, r3, r0    @ x[14] = x[14] ^ x[29]
vmov s14, r3    @ spill x[14] from r3
vmov r3, s3    @ load x[3] from s3
eor r0, r0, r3    @ x[29] = x[29] ^ x[3]
eor r0, r0, r7    @ x[29] = x[29] ^ x[13]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s8    @ load x[8] from s8
eor r2, r2, r4    @ x[5] = x[5] ^ x[8]
vmov s6, r1    @ spill x[6] from r1
vmov r1, s1    @ load x[1] from s1
eor r1, r1, r9    @ x[1] = x[1] ^ x[25]
eor r9, r9, r8    @ x[25] = x[25] ^ x[17]
vmov s1, r1    @ spill x[1] from r1
vmov r1, s7    @ load x[7] from s7
eor r9, r9, r1    @ x[25] = x[25] ^ x[7]
vmov s28, r5    @ spill x[28] from r5
vmov r5, s0    @ load x[0] from s0
vmov s29, r0    @ spill x[29] from r0
vmov r0, s10    @ load x[10] from s10
eor r5, r5, r0    @ x[0] = x[0] ^ x[10]
eor r3, r3, r10    @ x[3] = x[3] ^ x[23]
eor r10, r10, r7    @ x[23] = x[23] ^ x[13]
vmov s0, r5    @ spill x[0] from r5
vmov r5, s16    @ load x[16] from s16
eor r5, r5, r4    @ x[16] = x[16] ^ x[8]
vmov s30, r6    @ spill x[30] from r6
vmov r6, s9    @ load x[9] from s9
eor r5, r5, r6    @ x[16] = x[16] ^ x[9]
eor r6, r6, r7    @ x[9] = x[9] ^ x[13]
eor r3, r3, r0    @ x[3] = x[3] ^ x[10]
eor r7, r7, r8    @ x[13] = x[13] ^ x[17]
eor r2, r2, r0    @ x[5] = x[5] ^ x[10]
eor r4, r4, r1    @ x[8] = x[8] ^ x[7]
vmov s5, r2    @ spill x[5] from r2
vmov r2, s15    @ load x[15] from s15
eor r10, r10, r2    @ x[23] = x[23] ^ x[15]
vmov s7, r1    @ spill x[7] from r1
vmov r1, s20    @ load x[20] from s20
eor r2, r2, r1    @ x[15] = x[15] ^ x[20]
eor r11, r11, r7    @ x[11] = x[11] ^ x[13]
eor r10, r10, r9    @ x[23] = x[23] ^ x[25]
vmov s25, r9    @ spill x[25] from r9
vmov r9, s19    @ load x[19] from s19
eor r6, r6, r9    @ x[9] = x[9] ^ x[19]
vmov s11, r11    @ spill x[11] from r11
vmov r11, s12    @ load x[12] from s12
vmov s9, r6    @ spill x[9] from r6
vmov r6, s30    @ load x[30] from s30
eor r11, r11, r6    @ x[12] = x[12] ^ x[30]
eor r6, r6, r7    @ x[30] = x[30] ^ x[13]
eor r7, r7, r0    @ x[13] = x[13] ^ x[10]
vmov s10, r0    @ spill x[10] from r0
vmov r0, s29    @ load x[29] from s29
eor r6, r6, r0    @ x[30] = x[30] ^ x[29]
vmov s30, r6    @ spill x[30] from r6
vmov r6, s28    @ load x[28] from s28
eor r6, r6, r2    @ x[28] = x[28] ^ x[15]
eor r3, r3, r4    @ x[3] = x[3] ^ x[8]
vmov s3, r3    @ spill x[3] from r3
vmov r3, s0    @ load x[0] from s0
eor r4, r4, r3    @ x[8] = x[8] ^ x[0]
vmov s8, r4    @ spill x[8] from r4
vmov r4, s27    @ load x[27] from s27
eor r10, r10, r4    @ x[23] = x[23] ^ x[27]
vmov s28, r6    @ spill x[28] from r6
vmov r6, s6    @ load x[6] from s6
eor r6, r6, r2    @ x[6] = x[6] ^ x[15]
eor r0, r0, r1    @ x[29] = x[29] ^ x[20]
vmov s29, r0    @ spill x[29] from r0
vmov r0, s21    @ load x[21] from s21
eor r0, r0, r1    @ x[21] = x[21] ^ x[20]
eor r7, r7, r0    @ x[13] = x[13] ^ x[21]
vmov s20, r1    @ spill x[20] from r1
vmov r1, s26    @ load x[26] from s26
vmov s27, r4    @ spill x[27] from r4
vmov r4, s24    @ load x[24] from s24
eor r1, r1, r4    @ x[26] = x[26] ^ x[24]
eor r11, r11, r4    @ x[12] = x[12] ^ x[24]
eor r5, r5, r7    @ x[16] = x[16] ^ x[13]
vmov s13, r7    @ spill x[13] from r7
vmov r7, s4    @ load x[4] from s4
eor r7, r7, r11    @ x[4] = x[4] ^ x[12]
vmov s16, r5    @ spill x[16] from r5
vmov r5, s31    @ load x[31] from s31
eor r2, r2, r5    @ x[15] = x[15] ^ x[31]
eor r9, r9, r10    @ x[19] = x[19] ^ x[23]
vmov s23, r10    @ spill x[23] from r10
vmov r10, s1    @ load x[1] from s1
eor r3, r3, r10    @ x[0] = x[0] ^ x[1]
eor r8, r8, r5    @ x[17] = x[17] ^ x[31]
vmov s17, r8    @ spill x[17] from r8
vmov r8, s10    @ load x[10] from s10
eor r8, r8, r11    @ x[10] = x[10] ^ x[12]
vstr.32  s17, [r14, #0]    @ y[0] = x[17] from s17
vstr.32  s9, [r14, #4]    @ y[1] = x[9] from s9
str  r5, [r14, #8]    @ y[2] = x[31] from r5
vstr.32  s8, [r14, #12]    @ y[3] = x[8] from s8
vstr.32  s18, [r14, #16]    @ y[4] = x[18] from s18
str  r8, [r14, #20]    @ y[5] = x[10] from r8
vstr.32  s30, [r14, #24]    @ y[6] = x[30] from s30
str  r10, [r14, #28]    @ y[7] = x[1] from r10
str  r6, [r14, #32]    @ y[8] = x[6] from r6
vstr.32  s16, [r14, #36]    @ y[9] = x[16] from s16
vstr.32  s3, [r14, #40]    @ y[10] = x[3] from s3
vstr.32  s7, [r14, #44]    @ y[11] = x[7] from s7
vstr.32  s2, [r14, #48]    @ y[12] = x[2] from s2
vstr.32  s11, [r14, #52]    @ y[13] = x[11] from s11
str  r2, [r14, #56]    @ y[14] = x[15] from r2
vstr.32  s27, [r14, #60]    @ y[15] = x[27] from s27
str  r3, [r14, #64]    @ y[16] = x[0] from r3
str  r1, [r14, #68]    @ y[17] = x[26] from r1
str  r11, [r14, #72]    @ y[18] = x[12] from r11
vstr.32  s23, [r14, #76]    @ y[19] = x[23] from s23
vstr.32  s28, [r14, #80]    @ y[20] = x[28] from s28
vstr.32  s20, [r14, #84]    @ y[21] = x[20] from s20
vstr.32  s25, [r14, #88]    @ y[22] = x[25] from s25
str  r7, [r14, #92]    @ y[23] = x[4] from r7
vstr.32  s29, [r14, #96]    @ y[24] = x[29] from s29
vstr.32  s14, [r14, #100]    @ y[25] = x[14] from s14
str  r9, [r14, #104]    @ y[26] = x[19] from r9
vstr.32  s22, [r14, #108]    @ y[27] = x[22] from s22
vstr.32  s5, [r14, #112]    @ y[28] = x[5] from s5
vstr.32  s13, [r14, #116]    @ y[29] = x[13] from s13
str  r0, [r14, #120]    @ y[30] = x[21] from r0
str  r4, [r14, #124]    @ y[31] = x[24] from r4
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_v28, .-gft_mul_v28
