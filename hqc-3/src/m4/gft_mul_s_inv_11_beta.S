.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global gft_mul_s_inv_11_beta
.type gft_mul_s_inv_11_beta, %function
.align 2
gft_mul_s_inv_11_beta:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldm r12, {s0-s31}
vmov r0, s23    @ load x[23] from s23
vmov r1, s2    @ load x[2] from s2
eor r0, r0, r1    @ x[23] = x[23] ^ x[2]
vmov r2, s31    @ load x[31] from s31
vmov r3, s3    @ load x[3] from s3
eor r2, r2, r3    @ x[31] = x[31] ^ x[3]
vmov r4, s19    @ load x[19] from s19
eor r4, r4, r0    @ x[19] = x[19] ^ x[23]
vmov r5, s26    @ load x[26] from s26
eor r4, r4, r5    @ x[19] = x[19] ^ x[26]
vmov r6, s6    @ load x[6] from s6
vmov r7, s8    @ load x[8] from s8
eor r6, r6, r7    @ x[6] = x[6] ^ x[8]
vmov r8, s1    @ load x[1] from s1
eor r8, r8, r4    @ x[1] = x[1] ^ x[19]
vmov r9, s30    @ load x[30] from s30
eor r4, r4, r9    @ x[19] = x[19] ^ x[30]
vmov r10, s5    @ load x[5] from s5
vmov r11, s4    @ load x[4] from s4
eor r10, r10, r11    @ x[5] = x[5] ^ x[4]
vmov s1, r8    @ spill x[1] from r8
vmov r8, s16    @ load x[16] from s16
eor r8, r8, r10    @ x[16] = x[16] ^ x[5]
vmov s16, r8    @ spill x[16] from r8
vmov r8, s24    @ load x[24] from s24
vmov s30, r9    @ spill x[30] from r9
vmov r9, s13    @ load x[13] from s13
eor r8, r8, r9    @ x[24] = x[24] ^ x[13]
eor r3, r3, r6    @ x[3] = x[3] ^ x[6]
eor r5, r5, r3    @ x[26] = x[26] ^ x[3]
eor r7, r7, r2    @ x[8] = x[8] ^ x[31]
vmov s8, r7    @ spill x[8] from r7
vmov r7, s15    @ load x[15] from s15
eor r2, r2, r7    @ x[31] = x[31] ^ x[15]
vmov s6, r6    @ spill x[6] from r6
vmov r6, s14    @ load x[14] from s14
eor r0, r0, r6    @ x[23] = x[23] ^ x[14]
vmov s13, r9    @ spill x[13] from r9
vmov r9, s20    @ load x[20] from s20
eor r11, r11, r9    @ x[4] = x[4] ^ x[20]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s21    @ load x[21] from s21
eor r0, r0, r7    @ x[23] = x[23] ^ x[21]
vmov s3, r3    @ spill x[3] from r3
vmov r3, s18    @ load x[18] from s18
eor r7, r7, r3    @ x[21] = x[21] ^ x[18]
vmov s21, r7    @ spill x[21] from r7
vmov r7, s22    @ load x[22] from s22
vmov s18, r3    @ spill x[18] from r3
vmov r3, s25    @ load x[25] from s25
eor r7, r7, r3    @ x[22] = x[22] ^ x[25]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s12    @ load x[12] from s12
eor r3, r3, r4    @ x[12] = x[12] ^ x[19]
eor r1, r1, r0    @ x[2] = x[2] ^ x[23]
vmov s19, r4    @ spill x[19] from r4
vmov r4, s27    @ load x[27] from s27
eor r5, r5, r4    @ x[26] = x[26] ^ x[27]
vmov s26, r5    @ spill x[26] from r5
vmov r5, s10    @ load x[10] from s10
eor r5, r5, r7    @ x[10] = x[10] ^ x[22]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s30    @ load x[30] from s30
eor r10, r10, r7    @ x[5] = x[5] ^ x[30]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s0    @ load x[0] from s0
eor r2, r2, r4    @ x[31] = x[31] ^ x[0]
eor r10, r10, r6    @ x[5] = x[5] ^ x[14]
eor r9, r9, r3    @ x[20] = x[20] ^ x[12]
vmov s12, r3    @ spill x[12] from r3
vmov r3, s29    @ load x[29] from s29
vmov s14, r6    @ spill x[14] from r6
vmov r6, s7    @ load x[7] from s7
eor r3, r3, r6    @ x[29] = x[29] ^ x[7]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s28    @ load x[28] from s28
eor r6, r6, r2    @ x[28] = x[28] ^ x[31]
eor r7, r7, r0    @ x[30] = x[30] ^ x[23]
eor r2, r2, r11    @ x[31] = x[31] ^ x[4]
vmov s28, r6    @ spill x[28] from r6
vmov r6, s17    @ load x[17] from s17
eor r2, r2, r6    @ x[31] = x[31] ^ x[17]
vmov s31, r2    @ spill x[31] from r2
vmov r2, s11    @ load x[11] from s11
eor r2, r2, r10    @ x[11] = x[11] ^ x[5]
vmov s30, r7    @ spill x[30] from r7
vmov r7, s3    @ load x[3] from s3
eor r7, r7, r8    @ x[3] = x[3] ^ x[24]
eor r10, r10, r5    @ x[5] = x[5] ^ x[10]
eor r9, r9, r1    @ x[20] = x[20] ^ x[2]
eor r6, r6, r9    @ x[17] = x[17] ^ x[20]
vmov s3, r7    @ spill x[3] from r7
vmov r7, s9    @ load x[9] from s9
eor r5, r5, r7    @ x[10] = x[10] ^ x[9]
vmov s11, r2    @ spill x[11] from r2
vmov r2, s27    @ load x[27] from s27
eor r7, r7, r2    @ x[9] = x[9] ^ x[27]
eor r5, r5, r0    @ x[10] = x[10] ^ x[23]
vmov s9, r7    @ spill x[9] from r7
vmov r7, s15    @ load x[15] from s15
eor r7, r7, r6    @ x[15] = x[15] ^ x[17]
vmov s20, r9    @ spill x[20] from r9
vmov r9, s19    @ load x[19] from s19
eor r8, r8, r9    @ x[24] = x[24] ^ x[19]
vmov s10, r5    @ spill x[10] from r5
vmov r5, s26    @ load x[26] from s26
eor r9, r9, r5    @ x[19] = x[19] ^ x[26]
eor r2, r2, r3    @ x[27] = x[27] ^ x[29]
eor r2, r2, r4    @ x[27] = x[27] ^ x[0]
eor r0, r0, r10    @ x[23] = x[23] ^ x[5]
vmov s19, r9    @ spill x[19] from r9
vmov r9, s13    @ load x[13] from s13
eor r4, r4, r9    @ x[0] = x[0] ^ x[13]
vmov s23, r0    @ spill x[23] from r0
vmov r0, s18    @ load x[18] from s18
vmov s5, r10    @ spill x[5] from r10
vmov r10, s7    @ load x[7] from s7
eor r0, r0, r10    @ x[18] = x[18] ^ x[7]
eor r1, r1, r10    @ x[2] = x[2] ^ x[7]
eor r8, r8, r1    @ x[24] = x[24] ^ x[2]
eor r1, r1, r6    @ x[2] = x[2] ^ x[17]
eor r10, r10, r11    @ x[7] = x[7] ^ x[4]
vmov s7, r10    @ spill x[7] from r10
vmov r10, s14    @ load x[14] from s14
eor r11, r11, r10    @ x[4] = x[4] ^ x[14]
eor r7, r7, r11    @ x[15] = x[15] ^ x[4]
eor r10, r10, r9    @ x[14] = x[14] ^ x[13]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s10    @ load x[10] from s10
eor r9, r9, r7    @ x[13] = x[13] ^ x[10]
vmov s13, r9    @ spill x[13] from r9
vmov r9, s21    @ load x[21] from s21
eor r7, r7, r9    @ x[10] = x[10] ^ x[21]
eor r9, r9, r5    @ x[21] = x[21] ^ x[26]
eor r2, r2, r6    @ x[27] = x[27] ^ x[17]
vmov s17, r6    @ spill x[17] from r6
vmov r6, s20    @ load x[20] from s20
eor r3, r3, r6    @ x[29] = x[29] ^ x[20]
vmov s29, r3    @ spill x[29] from r3
vmov r3, s11    @ load x[11] from s11
eor r11, r11, r3    @ x[4] = x[4] ^ x[11]
vmov s11, r3    @ spill x[11] from r3
vmov r3, s9    @ load x[9] from s9
eor r6, r6, r3    @ x[20] = x[20] ^ x[9]
vmov s10, r7    @ spill x[10] from r7
vmov r7, s22    @ load x[22] from s22
eor r7, r7, r1    @ x[22] = x[22] ^ x[2]
eor r11, r11, r4    @ x[4] = x[4] ^ x[0]
vmov s4, r11    @ spill x[4] from r11
vmov r11, s30    @ load x[30] from s30
eor r11, r11, r7    @ x[30] = x[30] ^ x[22]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s6    @ load x[6] from s6
eor r7, r7, r2    @ x[6] = x[6] ^ x[27]
eor r10, r10, r7    @ x[14] = x[14] ^ x[6]
eor r3, r3, r0    @ x[9] = x[9] ^ x[18]
eor r5, r5, r6    @ x[26] = x[26] ^ x[20]
vmov s6, r7    @ spill x[6] from r7
vmov r7, s5    @ load x[5] from s5
eor r8, r8, r7    @ x[24] = x[24] ^ x[5]
eor r6, r6, r7    @ x[20] = x[20] ^ x[5]
eor r4, r4, r9    @ x[0] = x[0] ^ x[21]
vmov s0, r4    @ spill x[0] from r4
vmov r4, s23    @ load x[23] from s23
eor r9, r9, r4    @ x[21] = x[21] ^ x[23]
vmov s18, r0    @ spill x[18] from r0
vmov r0, s25    @ load x[25] from s25
eor r4, r4, r0    @ x[23] = x[23] ^ x[25]
eor r6, r6, r2    @ x[20] = x[20] ^ x[27]
vmov s20, r6    @ spill x[20] from r6
vmov r6, s7    @ load x[7] from s7
eor r7, r7, r6    @ x[5] = x[5] ^ x[7]
vmov s5, r7    @ spill x[5] from r7
vmov r7, s1    @ load x[1] from s1
eor r9, r9, r7    @ x[21] = x[21] ^ x[1]
eor r9, r9, r3    @ x[21] = x[21] ^ x[9]
eor r6, r6, r11    @ x[7] = x[7] ^ x[30]
vmov s24, r8    @ spill x[24] from r8
vmov r8, s12    @ load x[12] from s12
vmov s26, r5    @ spill x[26] from r5
vmov r5, s8    @ load x[8] from s8
eor r8, r8, r5    @ x[12] = x[12] ^ x[8]
vmov s27, r2    @ spill x[27] from r2
vmov r2, s3    @ load x[3] from s3
eor r7, r7, r2    @ x[1] = x[1] ^ x[3]
eor r11, r11, r3    @ x[30] = x[30] ^ x[9]
eor r9, r9, r1    @ x[21] = x[21] ^ x[2]
eor r6, r6, r2    @ x[7] = x[7] ^ x[3]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s10    @ load x[10] from s10
eor r2, r2, r11    @ x[3] = x[3] ^ x[10]
vmov s21, r9    @ spill x[21] from r9
vmov r9, s4    @ load x[4] from s4
eor r2, r2, r9    @ x[3] = x[3] ^ x[4]
vmov s3, r2    @ spill x[3] from r2
vmov r2, s19    @ load x[19] from s19
eor r9, r9, r2    @ x[4] = x[4] ^ x[19]
vmov s1, r7    @ spill x[1] from r7
vmov r7, s28    @ load x[28] from s28
eor r2, r2, r7    @ x[19] = x[19] ^ x[28]
eor r8, r8, r2    @ x[12] = x[12] ^ x[19]
eor r2, r2, r3    @ x[19] = x[19] ^ x[9]
vmov s12, r8    @ spill x[12] from r8
vmov r8, s11    @ load x[11] from s11
eor r3, r3, r8    @ x[9] = x[9] ^ x[11]
vmov s9, r3    @ spill x[9] from r3
vmov r3, s17    @ load x[17] from s17
eor r9, r9, r3    @ x[4] = x[4] ^ x[17]
eor r3, r3, r0    @ x[17] = x[17] ^ x[25]
vmov s17, r3    @ spill x[17] from r3
vmov r3, s29    @ load x[29] from s29
eor r0, r0, r3    @ x[25] = x[25] ^ x[29]
vmov s4, r9    @ spill x[4] from r9
vmov r9, s22    @ load x[22] from s22
eor r3, r3, r9    @ x[29] = x[29] ^ x[22]
vmov s25, r0    @ spill x[25] from r0
vmov r0, s13    @ load x[13] from s13
eor r0, r0, r8    @ x[13] = x[13] ^ x[11]
eor r1, r1, r8    @ x[2] = x[2] ^ x[11]
eor r8, r8, r10    @ x[11] = x[11] ^ x[14]
eor r10, r10, r9    @ x[14] = x[14] ^ x[22]
eor r9, r9, r0    @ x[22] = x[22] ^ x[13]
eor r9, r9, r1    @ x[22] = x[22] ^ x[2]
eor r11, r11, r4    @ x[10] = x[10] ^ x[23]
eor r1, r1, r5    @ x[2] = x[2] ^ x[8]
eor r4, r4, r5    @ x[23] = x[23] ^ x[8]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s16    @ load x[16] from s16
eor r5, r5, r11    @ x[8] = x[8] ^ x[16]
eor r8, r8, r6    @ x[11] = x[11] ^ x[7]
vmov s14, r10    @ spill x[14] from r10
vmov r10, s31    @ load x[31] from s31
eor r6, r6, r10    @ x[7] = x[7] ^ x[31]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s27    @ load x[27] from s27
eor r2, r2, r6    @ x[19] = x[19] ^ x[27]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s26    @ load x[26] from s26
eor r4, r4, r2    @ x[23] = x[23] ^ x[26]
eor r3, r3, r7    @ x[29] = x[29] ^ x[28]
vmov s29, r3    @ spill x[29] from r3
vmov r3, s15    @ load x[15] from s15
eor r9, r9, r3    @ x[22] = x[22] ^ x[15]
vmov s28, r7    @ spill x[28] from r7
vmov r7, s24    @ load x[24] from s24
eor r7, r7, r5    @ x[24] = x[24] ^ x[8]
vmov s26, r2    @ spill x[26] from r2
vmov r2, s1    @ load x[1] from s1
eor r6, r6, r2    @ x[27] = x[27] ^ x[1]
vmov s2, r1    @ spill x[2] from r1
vmov r1, s21    @ load x[21] from s21
eor r4, r4, r1    @ x[23] = x[23] ^ x[21]
vmov s21, r1    @ spill x[21] from r1
vmov r1, s25    @ load x[25] from s25
vmov s27, r6    @ spill x[27] from r6
vmov r6, s20    @ load x[20] from s20
eor r1, r1, r6    @ x[25] = x[25] ^ x[20]
eor r5, r5, r8    @ x[8] = x[8] ^ x[11]
eor r0, r0, r11    @ x[13] = x[13] ^ x[16]
eor r2, r2, r11    @ x[1] = x[1] ^ x[16]
eor r6, r6, r10    @ x[20] = x[20] ^ x[31]
vmov s20, r6    @ spill x[20] from r6
vmov r6, s18    @ load x[18] from s18
eor r10, r10, r6    @ x[31] = x[31] ^ x[18]
vmov s1, r2    @ spill x[1] from r2
vmov r2, s6    @ load x[6] from s6
eor r9, r9, r2    @ x[22] = x[22] ^ x[6]
eor r8, r8, r6    @ x[11] = x[11] ^ x[18]
eor r4, r4, r3    @ x[23] = x[23] ^ x[15]
eor r3, r3, r11    @ x[15] = x[15] ^ x[16]
eor r11, r11, r7    @ x[16] = x[16] ^ x[24]
vmov s23, r4    @ spill x[23] from r4
vmov r4, s27    @ load x[27] from s27
eor r2, r2, r4    @ x[6] = x[6] ^ x[27]
eor r7, r7, r4    @ x[24] = x[24] ^ x[27]
eor r7, r7, r10    @ x[24] = x[24] ^ x[31]
vmov s8, r5    @ spill x[8] from r5
vmov r5, s4    @ load x[4] from s4
eor r4, r4, r5    @ x[27] = x[27] ^ x[4]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s5    @ load x[5] from s5
vmov s22, r9    @ spill x[22] from r9
vmov r9, s2    @ load x[2] from s2
eor r4, r4, r9    @ x[5] = x[5] ^ x[2]
vmov s13, r0    @ spill x[13] from r0
vmov r0, s26    @ load x[26] from s26
vmov s2, r9    @ spill x[2] from r9
vmov r9, s12    @ load x[12] from s12
eor r0, r0, r9    @ x[26] = x[26] ^ x[12]
vmov s26, r0    @ spill x[26] from r0
vmov r0, s3    @ load x[3] from s3
eor r10, r10, r0    @ x[31] = x[31] ^ x[3]
eor r10, r10, r2    @ x[31] = x[31] ^ x[6]
vmov s15, r3    @ spill x[15] from r3
vmov r3, s14    @ load x[14] from s14
vmov s6, r2    @ spill x[6] from r2
vmov r2, s19    @ load x[19] from s19
eor r3, r3, r2    @ x[14] = x[14] ^ x[19]
vmov s11, r8    @ spill x[11] from r8
vmov r8, s17    @ load x[17] from s17
eor r2, r2, r8    @ x[19] = x[19] ^ x[17]
eor r8, r8, r6    @ x[17] = x[17] ^ x[18]
eor r6, r6, r9    @ x[18] = x[18] ^ x[12]
eor r5, r5, r10    @ x[4] = x[4] ^ x[31]
vmov s4, r5    @ spill x[4] from r5
vmov r5, s10    @ load x[10] from s10
eor r0, r0, r5    @ x[3] = x[3] ^ x[10]
eor r5, r5, r11    @ x[10] = x[10] ^ x[16]
vmov s16, r11    @ spill x[16] from r11
vmov r11, s28    @ load x[28] from s28
eor r0, r0, r11    @ x[3] = x[3] ^ x[28]
eor r11, r11, r6    @ x[28] = x[28] ^ x[18]
eor r6, r6, r3    @ x[18] = x[18] ^ x[14]
eor r11, r11, r4    @ x[28] = x[28] ^ x[5]
eor r10, r10, r7    @ x[31] = x[31] ^ x[24]
vmov s3, r0    @ spill x[3] from r0
vmov r0, s30    @ load x[30] from s30
eor r3, r3, r0    @ x[14] = x[14] ^ x[30]
eor r0, r0, r1    @ x[30] = x[30] ^ x[25]
vmov s30, r0    @ spill x[30] from r0
vmov r0, s21    @ load x[21] from s21
eor r10, r10, r0    @ x[31] = x[31] ^ x[21]
eor r0, r0, r1    @ x[21] = x[21] ^ x[25]
eor r9, r9, r0    @ x[12] = x[12] ^ x[21]
eor r1, r1, r4    @ x[25] = x[25] ^ x[5]
vmov s12, r9    @ spill x[12] from r9
vmov r9, s11    @ load x[11] from s11
eor r4, r4, r9    @ x[5] = x[5] ^ x[11]
eor r0, r0, r2    @ x[21] = x[21] ^ x[19]
eor r9, r9, r2    @ x[11] = x[11] ^ x[19]
eor r2, r2, r5    @ x[19] = x[19] ^ x[10]
eor r5, r5, r8    @ x[10] = x[10] ^ x[17]
vmov s31, r10    @ spill x[31] from r10
vmov r10, s6    @ load x[6] from s6
eor r9, r9, r10    @ x[11] = x[11] ^ x[6]
vmov s19, r2    @ spill x[19] from r2
vmov r2, s15    @ load x[15] from s15
eor r10, r10, r2    @ x[6] = x[6] ^ x[15]
vmov s6, r10    @ spill x[6] from r10
vmov r10, s9    @ load x[9] from s9
eor r2, r2, r10    @ x[15] = x[15] ^ x[9]
eor r10, r10, r6    @ x[9] = x[9] ^ x[18]
eor r7, r7, r11    @ x[24] = x[24] ^ x[28]
eor r6, r6, r7    @ x[18] = x[18] ^ x[24]
vmov s28, r11    @ spill x[28] from r11
vmov r11, s29    @ load x[29] from s29
eor r8, r8, r11    @ x[17] = x[17] ^ x[29]
eor r9, r9, r11    @ x[11] = x[11] ^ x[29]
eor r7, r7, r11    @ x[24] = x[24] ^ x[29]
vmov s17, r8    @ spill x[17] from r8
vmov r8, s2    @ load x[2] from s2
eor r10, r10, r8    @ x[9] = x[9] ^ x[2]
vmov s2, r8    @ spill x[2] from r8
vmov r8, s16    @ load x[16] from s16
eor r11, r11, r8    @ x[29] = x[29] ^ x[16]
eor r2, r2, r8    @ x[15] = x[15] ^ x[16]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s13    @ load x[13] from s13
eor r8, r8, r6    @ x[16] = x[16] ^ x[13]
vmov s16, r8    @ spill x[16] from r8
vmov r8, s22    @ load x[22] from s22
eor r7, r7, r8    @ x[24] = x[24] ^ x[22]
vmov s13, r6    @ spill x[13] from r6
vmov r6, s8    @ load x[8] from s8
eor r3, r3, r6    @ x[14] = x[14] ^ x[8]
eor r4, r4, r6    @ x[5] = x[5] ^ x[8]
vmov s22, r8    @ spill x[22] from r8
vmov r8, s7    @ load x[7] from s7
vmov s8, r6    @ spill x[8] from r6
vmov r6, s1    @ load x[1] from s1
eor r8, r8, r6    @ x[7] = x[7] ^ x[1]
vmov s7, r8    @ spill x[7] from r8
vmov r8, s0    @ load x[0] from s0
vmov s11, r9    @ spill x[11] from r9
vmov r9, s20    @ load x[20] from s20
eor r8, r8, r9    @ x[0] = x[0] ^ x[20]
vmov s0, r8    @ spill x[0] from r8
vmov r8, s23    @ load x[23] from s23
eor r8, r8, r11    @ x[23] = x[23] ^ x[29]
eor r3, r3, r0    @ x[14] = x[14] ^ x[21]
vmov s14, r3    @ spill x[14] from r3
vmov r3, s4    @ load x[4] from s4
eor r3, r3, r1    @ x[4] = x[4] ^ x[25]
eor r1, r1, r2    @ x[25] = x[25] ^ x[15]
eor r2, r2, r4    @ x[15] = x[15] ^ x[5]
eor r6, r6, r5    @ x[1] = x[1] ^ x[10]
vmov s5, r4    @ spill x[5] from r4
vmov r4, s19    @ load x[19] from s19
eor r4, r4, r3    @ x[19] = x[19] ^ x[4]
vmov s10, r5    @ spill x[10] from r5
vmov r5, s26    @ load x[26] from s26
eor r3, r3, r5    @ x[4] = x[4] ^ x[26]
vmov s21, r0    @ spill x[21] from r0
vmov r0, s27    @ load x[27] from s27
vmov s15, r2    @ spill x[15] from r2
vmov r2, s31    @ load x[31] from s31
eor r0, r0, r2    @ x[27] = x[27] ^ x[31]
eor r7, r7, r1    @ x[24] = x[24] ^ x[25]
vmov s29, r11    @ spill x[29] from r11
vmov r11, s12    @ load x[12] from s12
eor r11, r11, r10    @ x[12] = x[12] ^ x[9]
vmov s9, r10    @ spill x[9] from r10
vmov r10, s11    @ load x[11] from s11
eor r1, r1, r10    @ x[25] = x[25] ^ x[11]
vmov s25, r1    @ spill x[25] from r1
vmov r1, s18    @ load x[18] from s18
eor r1, r1, r8    @ x[18] = x[18] ^ x[23]
vmov s1, r6    @ spill x[1] from r6
vmov r6, s30    @ load x[30] from s30
eor r6, r6, r7    @ x[30] = x[30] ^ x[24]
vmov s24, r7    @ spill x[24] from r7
vmov r7, s2    @ load x[2] from s2
eor r7, r7, r1    @ x[2] = x[2] ^ x[18]
vstr.32  s16, [r14, #0]    @ y[0] = x[16] from s16
vstr.32  s28, [r14, #4]    @ y[1] = x[28] from s28
vstr.32  s8, [r14, #8]    @ y[2] = x[8] from s8
str  r6, [r14, #12]    @ y[3] = x[30] from r6
vstr.32  s6, [r14, #16]    @ y[4] = x[6] from s6
vstr.32  s7, [r14, #20]    @ y[5] = x[7] from s7
vstr.32  s17, [r14, #24]    @ y[6] = x[17] from s17
str  r11, [r14, #28]    @ y[7] = x[12] from r11
vstr.32  s5, [r14, #32]    @ y[8] = x[5] from s5
str  r0, [r14, #36]    @ y[9] = x[27] from r0
vstr.32  s24, [r14, #40]    @ y[10] = x[24] from s24
str  r1, [r14, #44]    @ y[11] = x[18] from r1
vstr.32  s13, [r14, #48]    @ y[12] = x[13] from s13
vstr.32  s10, [r14, #52]    @ y[13] = x[10] from s10
vstr.32  s21, [r14, #56]    @ y[14] = x[21] from s21
str  r5, [r14, #60]    @ y[15] = x[26] from r5
str  r7, [r14, #64]    @ y[16] = x[2] from r7
str  r10, [r14, #68]    @ y[17] = x[11] from r10
vstr.32  s0, [r14, #72]    @ y[18] = x[0] from s0
str  r4, [r14, #76]    @ y[19] = x[19] from r4
vstr.32  s25, [r14, #80]    @ y[20] = x[25] from s25
vstr.32  s29, [r14, #84]    @ y[21] = x[29] from s29
vstr.32  s1, [r14, #88]    @ y[22] = x[1] from s1
vstr.32  s14, [r14, #92]    @ y[23] = x[14] from s14
vstr.32  s3, [r14, #96]    @ y[24] = x[3] from s3
str  r9, [r14, #100]    @ y[25] = x[20] from r9
str  r3, [r14, #104]    @ y[26] = x[4] from r3
vstr.32  s22, [r14, #108]    @ y[27] = x[22] from s22
vstr.32  s15, [r14, #112]    @ y[28] = x[15] from s15
vstr.32  s9, [r14, #116]    @ y[29] = x[9] from s9
str  r2, [r14, #120]    @ y[30] = x[31] from r2
str  r8, [r14, #124]    @ y[31] = x[23] from r8
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   gft_mul_s_inv_11_beta, .-gft_mul_s_inv_11_beta