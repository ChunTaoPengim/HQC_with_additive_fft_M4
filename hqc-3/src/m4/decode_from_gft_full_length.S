.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global decode_from_gft_full_length
.type decode_from_gft_full_length, %function
.align 2
decode_from_gft_full_length:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldm r12, {s0-s31}
vmov r0, s30    @ load x[30] from s30
vmov r1, s26    @ load x[26] from s26
eor r0, r0, r1    @ x[30] = x[30] ^ x[26]
vmov r2, s12    @ load x[12] from s12
eor r2, r2, r0    @ x[12] = x[12] ^ x[30]
vmov r3, s1    @ load x[1] from s1
vmov r4, s31    @ load x[31] from s31
eor r3, r3, r4    @ x[1] = x[1] ^ x[31]
vmov r5, s17    @ load x[17] from s17
vmov r6, s23    @ load x[23] from s23
eor r5, r5, r6    @ x[17] = x[17] ^ x[23]
vmov r7, s7    @ load x[7] from s7
vmov r8, s9    @ load x[9] from s9
eor r7, r7, r8    @ x[7] = x[7] ^ x[9]
vmov r9, s25    @ load x[25] from s25
eor r9, r9, r7    @ x[25] = x[25] ^ x[7]
vmov r10, s19    @ load x[19] from s19
vmov r11, s15    @ load x[15] from s15
eor r10, r10, r11    @ x[19] = x[19] ^ x[15]
vmov s31, r4    @ spill x[31] from r4
vmov r4, s27    @ load x[27] from s27
vmov s7, r7    @ spill x[7] from r7
vmov r7, s22    @ load x[22] from s22
eor r4, r4, r7    @ x[27] = x[27] ^ x[22]
eor r7, r7, r2    @ x[22] = x[22] ^ x[12]
vmov s25, r9    @ spill x[25] from r9
vmov r9, s14    @ load x[14] from s14
eor r9, r9, r7    @ x[14] = x[14] ^ x[22]
vmov s12, r2    @ spill x[12] from r2
vmov r2, s11    @ load x[11] from s11
eor r2, r2, r5    @ x[11] = x[11] ^ x[17]
vmov s30, r0    @ spill x[30] from r0
vmov r0, s18    @ load x[18] from s18
vmov s19, r10    @ spill x[19] from r10
vmov r10, s21    @ load x[21] from s21
eor r0, r0, r10    @ x[18] = x[18] ^ x[21]
eor r11, r11, r5    @ x[15] = x[15] ^ x[17]
vmov s15, r11    @ spill x[15] from r11
vmov r11, s10    @ load x[10] from s10
eor r1, r1, r11    @ x[26] = x[26] ^ x[10]
eor r10, r10, r7    @ x[21] = x[21] ^ x[22]
eor r1, r1, r0    @ x[26] = x[26] ^ x[18]
eor r4, r4, r2    @ x[27] = x[27] ^ x[11]
vmov s18, r0    @ spill x[18] from r0
vmov r0, s2    @ load x[2] from s2
eor r0, r0, r9    @ x[2] = x[2] ^ x[14]
eor r5, r5, r3    @ x[17] = x[17] ^ x[1]
eor r5, r5, r9    @ x[17] = x[17] ^ x[14]
vmov s17, r5    @ spill x[17] from r5
vmov r5, s29    @ load x[29] from s29
vmov s10, r11    @ spill x[10] from r11
vmov r11, s24    @ load x[24] from s24
eor r5, r5, r11    @ x[29] = x[29] ^ x[24]
eor r11, r11, r10    @ x[24] = x[24] ^ x[21]
vmov s29, r5    @ spill x[29] from r5
vmov r5, s5    @ load x[5] from s5
eor r5, r5, r3    @ x[5] = x[5] ^ x[1]
vmov s24, r11    @ spill x[24] from r11
vmov r11, s8    @ load x[8] from s8
eor r11, r11, r5    @ x[8] = x[8] ^ x[5]
vmov s11, r2    @ spill x[11] from r2
vmov r2, s3    @ load x[3] from s3
eor r8, r8, r2    @ x[9] = x[9] ^ x[3]
eor r2, r2, r11    @ x[3] = x[3] ^ x[8]
vmov s9, r8    @ spill x[9] from r8
vmov r8, s28    @ load x[28] from s28
eor r6, r6, r8    @ x[23] = x[23] ^ x[28]
eor r2, r2, r3    @ x[3] = x[3] ^ x[1]
vmov s23, r6    @ spill x[23] from r6
vmov r6, s19    @ load x[19] from s19
eor r3, r3, r6    @ x[1] = x[1] ^ x[19]
eor r6, r6, r9    @ x[19] = x[19] ^ x[14]
eor r6, r6, r1    @ x[19] = x[19] ^ x[26]
eor r1, r1, r11    @ x[26] = x[26] ^ x[8]
vmov s19, r6    @ spill x[19] from r6
vmov r6, s13    @ load x[13] from s13
eor r8, r8, r6    @ x[28] = x[28] ^ x[13]
vmov s1, r3    @ spill x[1] from r3
vmov r3, s30    @ load x[30] from s30
eor r9, r9, r3    @ x[14] = x[14] ^ x[30]
eor r9, r9, r6    @ x[14] = x[14] ^ x[13]
eor r8, r8, r0    @ x[28] = x[28] ^ x[2]
vmov s28, r8    @ spill x[28] from r8
vmov r8, s16    @ load x[16] from s16
eor r3, r3, r8    @ x[30] = x[30] ^ x[16]
vmov s30, r3    @ spill x[30] from r3
vmov r3, s12    @ load x[12] from s12
eor r6, r6, r3    @ x[13] = x[13] ^ x[12]
eor r11, r11, r6    @ x[8] = x[8] ^ x[13]
eor r1, r1, r8    @ x[26] = x[26] ^ x[16]
vmov s26, r1    @ spill x[26] from r1
vmov r1, s4    @ load x[4] from s4
eor r6, r6, r1    @ x[13] = x[13] ^ x[4]
eor r10, r10, r3    @ x[21] = x[21] ^ x[12]
eor r8, r8, r3    @ x[16] = x[16] ^ x[12]
vmov s21, r10    @ spill x[21] from r10
vmov r10, s0    @ load x[0] from s0
vmov s8, r11    @ spill x[8] from r11
vmov r11, s15    @ load x[15] from s15
eor r10, r10, r11    @ x[0] = x[0] ^ x[15]
eor r5, r5, r4    @ x[5] = x[5] ^ x[27]
eor r8, r8, r11    @ x[16] = x[16] ^ x[15]
eor r11, r11, r1    @ x[15] = x[15] ^ x[4]
vmov s15, r11    @ spill x[15] from r11
vmov r11, s6    @ load x[6] from s6
eor r7, r7, r11    @ x[22] = x[22] ^ x[6]
vmov s16, r8    @ spill x[16] from r8
vmov r8, s25    @ load x[25] from s25
eor r0, r0, r8    @ x[2] = x[2] ^ x[25]
vmov s2, r0    @ spill x[2] from r0
vmov r0, s11    @ load x[11] from s11
eor r3, r3, r0    @ x[12] = x[12] ^ x[11]
vmov s5, r5    @ spill x[5] from r5
vmov r5, s10    @ load x[10] from s10
eor r4, r4, r5    @ x[27] = x[27] ^ x[10]
eor r9, r9, r2    @ x[14] = x[14] ^ x[3]
eor r2, r2, r5    @ x[3] = x[3] ^ x[10]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s20    @ load x[20] from s20
eor r5, r5, r4    @ x[10] = x[10] ^ x[20]
vmov s20, r4    @ spill x[20] from r4
vmov r4, s24    @ load x[24] from s24
vmov s10, r5    @ spill x[10] from r5
vmov r5, s1    @ load x[1] from s1
eor r4, r4, r5    @ x[24] = x[24] ^ x[1]
eor r0, r0, r5    @ x[11] = x[11] ^ x[1]
eor r10, r10, r2    @ x[0] = x[0] ^ x[3]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s17    @ load x[17] from s17
eor r4, r4, r2    @ x[17] = x[17] ^ x[3]
eor r5, r5, r8    @ x[1] = x[1] ^ x[25]
eor r1, r1, r8    @ x[4] = x[4] ^ x[25]
eor r8, r8, r2    @ x[25] = x[25] ^ x[3]
eor r2, r2, r11    @ x[3] = x[3] ^ x[6]
eor r11, r11, r6    @ x[6] = x[6] ^ x[13]
vmov s6, r11    @ spill x[6] from r11
vmov r11, s30    @ load x[30] from s30
eor r6, r6, r11    @ x[13] = x[13] ^ x[30]
vmov s13, r6    @ spill x[13] from r6
vmov r6, s29    @ load x[29] from s29
eor r11, r11, r6    @ x[30] = x[30] ^ x[29]
eor r3, r3, r6    @ x[12] = x[12] ^ x[29]
eor r9, r9, r7    @ x[14] = x[14] ^ x[22]
vmov s4, r1    @ spill x[4] from r1
vmov r1, s23    @ load x[23] from s23
eor r1, r1, r7    @ x[23] = x[23] ^ x[22]
vmov s14, r9    @ spill x[14] from r9
vmov r9, s8    @ load x[8] from s8
eor r2, r2, r9    @ x[3] = x[3] ^ x[8]
vmov s3, r2    @ spill x[3] from r2
vmov r2, s7    @ load x[7] from s7
eor r9, r9, r2    @ x[8] = x[8] ^ x[7]
vmov s8, r9    @ spill x[8] from r9
vmov r9, s28    @ load x[28] from s28
eor r6, r6, r9    @ x[29] = x[29] ^ x[28]
vmov s29, r6    @ spill x[29] from r6
vmov r6, s27    @ load x[27] from s27
eor r9, r9, r6    @ x[28] = x[28] ^ x[27]
eor r10, r10, r6    @ x[0] = x[0] ^ x[27]
vmov s0, r10    @ spill x[0] from r10
vmov r10, s5    @ load x[5] from s5
eor r10, r10, r2    @ x[5] = x[5] ^ x[7]
vmov s28, r9    @ spill x[28] from r9
vmov r9, s26    @ load x[26] from s26
eor r9, r9, r6    @ x[26] = x[26] ^ x[27]
vmov s26, r9    @ spill x[26] from r9
vmov r9, s10    @ load x[10] from s10
eor r2, r2, r9    @ x[7] = x[7] ^ x[10]
vmov s11, r0    @ spill x[11] from r0
vmov r0, s2    @ load x[2] from s2
eor r6, r6, r0    @ x[27] = x[27] ^ x[2]
vmov s27, r6    @ spill x[27] from r6
vmov r6, s21    @ load x[21] from s21
eor r9, r9, r6    @ x[10] = x[10] ^ x[21]
eor r1, r1, r0    @ x[23] = x[23] ^ x[2]
eor r6, r6, r8    @ x[21] = x[21] ^ x[25]
vmov s23, r1    @ spill x[23] from r1
vmov r1, s18    @ load x[18] from s18
eor r8, r8, r1    @ x[25] = x[25] ^ x[18]
vmov s21, r6    @ spill x[21] from r6
vmov r6, s16    @ load x[16] from s16
eor r6, r6, r0    @ x[16] = x[16] ^ x[2]
eor r11, r11, r0    @ x[30] = x[30] ^ x[2]
eor r1, r1, r10    @ x[18] = x[18] ^ x[5]
vmov s30, r11    @ spill x[30] from r11
vmov r11, s9    @ load x[9] from s9
eor r4, r4, r11    @ x[17] = x[17] ^ x[9]
eor r0, r0, r5    @ x[2] = x[2] ^ x[1]
eor r6, r6, r8    @ x[16] = x[16] ^ x[25]
eor r10, r10, r6    @ x[5] = x[5] ^ x[16]
eor r6, r6, r11    @ x[16] = x[16] ^ x[9]
eor r2, r2, r11    @ x[7] = x[7] ^ x[9]
eor r5, r5, r7    @ x[1] = x[1] ^ x[22]
vmov s16, r6    @ spill x[16] from r6
vmov r6, s20    @ load x[20] from s20
eor r6, r6, r1    @ x[20] = x[20] ^ x[18]
eor r3, r3, r6    @ x[12] = x[12] ^ x[20]
vmov s12, r3    @ spill x[12] from r3
vmov r3, s14    @ load x[14] from s14
vmov s1, r5    @ spill x[1] from r5
vmov r5, s11    @ load x[11] from s11
eor r3, r3, r5    @ x[14] = x[14] ^ x[11]
vmov s17, r4    @ spill x[17] from r4
vmov r4, s4    @ load x[4] from s4
eor r4, r4, r10    @ x[4] = x[4] ^ x[5]
eor r10, r10, r6    @ x[5] = x[5] ^ x[20]
vmov s25, r8    @ spill x[25] from r8
vmov r8, s31    @ load x[31] from s31
eor r9, r9, r8    @ x[10] = x[10] ^ x[31]
eor r7, r7, r9    @ x[22] = x[22] ^ x[10]
vmov s10, r9    @ spill x[10] from r9
vmov r9, s13    @ load x[13] from s13
eor r2, r2, r9    @ x[7] = x[7] ^ x[13]
vmov s7, r2    @ spill x[7] from r2
vmov r2, s28    @ load x[28] from s28
eor r0, r0, r2    @ x[2] = x[2] ^ x[28]
vmov s2, r0    @ spill x[2] from r0
vmov r0, s6    @ load x[6] from s6
vmov s28, r2    @ spill x[28] from r2
vmov r2, s21    @ load x[21] from s21
eor r0, r0, r2    @ x[6] = x[6] ^ x[21]
vmov s6, r0    @ spill x[6] from r0
vmov r0, s27    @ load x[27] from s27
eor r2, r2, r0    @ x[21] = x[21] ^ x[27]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s0    @ load x[0] from s0
eor r7, r7, r3    @ x[0] = x[0] ^ x[14]
eor r0, r0, r1    @ x[27] = x[27] ^ x[18]
vmov s31, r8    @ spill x[31] from r8
vmov r8, s8    @ load x[8] from s8
eor r1, r1, r8    @ x[18] = x[18] ^ x[8]
vmov s18, r1    @ spill x[18] from r1
vmov r1, s19    @ load x[19] from s19
eor r1, r1, r10    @ x[19] = x[19] ^ x[5]
eor r3, r3, r6    @ x[14] = x[14] ^ x[20]
eor r4, r4, r5    @ x[4] = x[4] ^ x[11]
vmov s19, r1    @ spill x[19] from r1
vmov r1, s29    @ load x[29] from s29
eor r1, r1, r11    @ x[29] = x[29] ^ x[9]
vmov s4, r4    @ spill x[4] from r4
vmov r4, s26    @ load x[26] from s26
eor r9, r9, r4    @ x[13] = x[13] ^ x[26]
eor r7, r7, r8    @ x[0] = x[0] ^ x[8]
vmov s0, r7    @ spill x[0] from r7
vmov r7, s15    @ load x[15] from s15
eor r8, r8, r7    @ x[8] = x[8] ^ x[15]
vmov s8, r8    @ spill x[8] from r8
vmov r8, s30    @ load x[30] from s30
vmov s13, r9    @ spill x[13] from r9
vmov r9, s25    @ load x[25] from s25
eor r8, r8, r9    @ x[30] = x[30] ^ x[25]
eor r0, r0, r5    @ x[27] = x[27] ^ x[11]
eor r5, r5, r2    @ x[11] = x[11] ^ x[21]
vmov s30, r8    @ spill x[30] from r8
vmov r8, s31    @ load x[31] from s31
eor r9, r9, r8    @ x[25] = x[25] ^ x[31]
vmov s27, r0    @ spill x[27] from r0
vmov r0, s22    @ load x[22] from s22
eor r11, r11, r0    @ x[9] = x[9] ^ x[22]
eor r0, r0, r7    @ x[22] = x[22] ^ x[15]
eor r8, r8, r7    @ x[31] = x[31] ^ x[15]
vmov s22, r0    @ spill x[22] from r0
vmov r0, s3    @ load x[3] from s3
eor r7, r7, r0    @ x[15] = x[15] ^ x[3]
vmov s9, r11    @ spill x[9] from r11
vmov r11, s17    @ load x[17] from s17
eor r7, r7, r11    @ x[15] = x[15] ^ x[17]
vmov s17, r11    @ spill x[17] from r11
vmov r11, s10    @ load x[10] from s10
eor r3, r3, r11    @ x[14] = x[14] ^ x[10]
eor r6, r6, r8    @ x[20] = x[20] ^ x[31]
eor r8, r8, r4    @ x[31] = x[31] ^ x[26]
eor r1, r1, r0    @ x[29] = x[29] ^ x[3]
vmov s20, r6    @ spill x[20] from r6
vmov r6, s24    @ load x[24] from s24
eor r0, r0, r6    @ x[3] = x[3] ^ x[24]
vmov s14, r3    @ spill x[14] from r3
vmov r3, s1    @ load x[1] from s1
eor r9, r9, r3    @ x[25] = x[25] ^ x[1]
eor r9, r9, r0    @ x[25] = x[25] ^ x[3]
eor r5, r5, r6    @ x[11] = x[11] ^ x[24]
vmov s29, r1    @ spill x[29] from r1
vmov r1, s28    @ load x[28] from s28
eor r6, r6, r1    @ x[24] = x[24] ^ x[28]
vmov s25, r9    @ spill x[25] from r9
vmov r9, s6    @ load x[6] from s6
eor r1, r1, r9    @ x[28] = x[28] ^ x[6]
eor r0, r0, r10    @ x[3] = x[3] ^ x[5]
vmov s31, r8    @ spill x[31] from r8
vmov r8, s2    @ load x[2] from s2
eor r10, r10, r8    @ x[5] = x[5] ^ x[2]
eor r8, r8, r4    @ x[2] = x[2] ^ x[26]
eor r4, r4, r9    @ x[26] = x[26] ^ x[6]
eor r4, r4, r3    @ x[26] = x[26] ^ x[1]
eor r9, r9, r2    @ x[6] = x[6] ^ x[21]
vmov s6, r9    @ spill x[6] from r9
vmov r9, s7    @ load x[7] from s7
eor r1, r1, r9    @ x[28] = x[28] ^ x[7]
vmov s28, r1    @ spill x[28] from r1
vmov r1, s23    @ load x[23] from s23
eor r9, r9, r1    @ x[7] = x[7] ^ x[23]
eor r2, r2, r7    @ x[21] = x[21] ^ x[15]
vmov s7, r9    @ spill x[7] from r9
vmov r9, s27    @ load x[27] from s27
vmov s21, r2    @ spill x[21] from r2
vmov r2, s12    @ load x[12] from s12
eor r9, r9, r2    @ x[27] = x[27] ^ x[12]
eor r7, r7, r1    @ x[15] = x[15] ^ x[23]
eor r1, r1, r3    @ x[23] = x[23] ^ x[1]
vmov s27, r9    @ spill x[27] from r9
vmov r9, s17    @ load x[17] from s17
eor r3, r3, r9    @ x[1] = x[1] ^ x[17]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s18    @ load x[18] from s18
eor r11, r11, r7    @ x[10] = x[10] ^ x[18]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s13    @ load x[13] from s13
eor r9, r9, r11    @ x[17] = x[17] ^ x[13]
vmov s17, r9    @ spill x[17] from r9
vmov r9, s9    @ load x[9] from s9
eor r9, r9, r3    @ x[9] = x[9] ^ x[1]
eor r3, r3, r2    @ x[1] = x[1] ^ x[12]
eor r0, r0, r4    @ x[3] = x[3] ^ x[26]
vmov s9, r9    @ spill x[9] from r9
vmov r9, s8    @ load x[8] from s8
eor r4, r4, r9    @ x[26] = x[26] ^ x[8]
vmov s26, r4    @ spill x[26] from r4
vmov r4, s16    @ load x[16] from s16
eor r10, r10, r4    @ x[5] = x[5] ^ x[16]
eor r0, r0, r11    @ x[3] = x[3] ^ x[13]
vmov s3, r0    @ spill x[3] from r0
vmov r0, s22    @ load x[22] from s22
eor r0, r0, r6    @ x[22] = x[22] ^ x[24]
eor r5, r5, r1    @ x[11] = x[11] ^ x[23]
vmov s11, r5    @ spill x[11] from r5
vmov r5, s31    @ load x[31] from s31
eor r0, r0, r5    @ x[22] = x[22] ^ x[31]
eor r7, r7, r1    @ x[18] = x[18] ^ x[23]
vmov s22, r0    @ spill x[22] from r0
vmov r0, s4    @ load x[4] from s4
eor r0, r0, r8    @ x[4] = x[4] ^ x[2]
eor r9, r9, r6    @ x[8] = x[8] ^ x[24]
vmov s2, r8    @ spill x[2] from r8
vmov r8, s15    @ load x[15] from s15
eor r11, r11, r8    @ x[13] = x[13] ^ x[15]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s25    @ load x[25] from s25
eor r6, r6, r11    @ x[24] = x[24] ^ x[25]
eor r8, r8, r5    @ x[15] = x[15] ^ x[31]
eor r1, r1, r10    @ x[23] = x[23] ^ x[5]
vmov s24, r6    @ spill x[24] from r6
vmov r6, s21    @ load x[21] from s21
eor r4, r4, r6    @ x[16] = x[16] ^ x[21]
vmov s16, r4    @ spill x[16] from r4
vmov r4, s17    @ load x[17] from s17
eor r11, r11, r4    @ x[25] = x[25] ^ x[17]
vmov s15, r8    @ spill x[15] from r8
vmov r8, s29    @ load x[29] from s29
eor r1, r1, r8    @ x[23] = x[23] ^ x[29]
vmov s23, r1    @ spill x[23] from r1
vmov r1, s14    @ load x[14] from s14
vmov s21, r6    @ spill x[21] from r6
vmov r6, s30    @ load x[30] from s30
eor r1, r1, r6    @ x[14] = x[14] ^ x[30]
vmov s29, r8    @ spill x[29] from r8
vmov r8, s6    @ load x[6] from s6
eor r5, r5, r8    @ x[31] = x[31] ^ x[6]
eor r4, r4, r2    @ x[17] = x[17] ^ x[12]
vmov s12, r2    @ spill x[12] from r2
vmov r2, s20    @ load x[20] from s20
eor r2, r2, r10    @ x[20] = x[20] ^ x[5]
vmov s5, r10    @ spill x[5] from r10
vmov r10, s10    @ load x[10] from s10
eor r11, r11, r10    @ x[25] = x[25] ^ x[10]
vmov s31, r5    @ spill x[31] from r5
vmov r5, s27    @ load x[27] from s27
eor r7, r7, r5    @ x[18] = x[18] ^ x[27]
vmov s18, r7    @ spill x[18] from r7
vmov r7, s19    @ load x[19] from s19
eor r7, r7, r0    @ x[19] = x[19] ^ x[4]
eor r1, r1, r9    @ x[14] = x[14] ^ x[8]
eor r3, r3, r11    @ x[1] = x[1] ^ x[25]
vmov s8, r9    @ spill x[8] from r9
vmov r9, s13    @ load x[13] from s13
eor r9, r9, r3    @ x[13] = x[13] ^ x[1]
vmov s1, r3    @ spill x[1] from r3
vmov r3, s2    @ load x[2] from s2
eor r3, r3, r11    @ x[2] = x[2] ^ x[25]
eor r8, r8, r9    @ x[6] = x[6] ^ x[13]
eor r10, r10, r9    @ x[10] = x[10] ^ x[13]
vmov s4, r0    @ spill x[4] from r0
vmov r0, s22    @ load x[22] from s22
eor r8, r8, r0    @ x[6] = x[6] ^ x[22]
vstr.32  s0, [r14, #0]    @ y[0] = x[0] from s0
str  r11, [r14, #4]    @ y[1] = x[25] from r11
str  r1, [r14, #8]    @ y[2] = x[14] from r1
vstr.32  s9, [r14, #12]    @ y[3] = x[9] from s9
vstr.32  s5, [r14, #16]    @ y[4] = x[5] from s5
str  r9, [r14, #20]    @ y[5] = x[13] from r9
vstr.32  s11, [r14, #24]    @ y[6] = x[11] from s11
vstr.32  s28, [r14, #28]    @ y[7] = x[28] from s28
str  r7, [r14, #32]    @ y[8] = x[19] from r7
str  r5, [r14, #36]    @ y[9] = x[27] from r5
vstr.32  s1, [r14, #40]    @ y[10] = x[1] from s1
vstr.32  s7, [r14, #44]    @ y[11] = x[7] from s7
str  r4, [r14, #48]    @ y[12] = x[17] from r4
str  r0, [r14, #52]    @ y[13] = x[22] from r0
vstr.32  s24, [r14, #56]    @ y[14] = x[24] from s24
vstr.32  s4, [r14, #60]    @ y[15] = x[4] from s4
vstr.32  s31, [r14, #64]    @ y[16] = x[31] from s31
vstr.32  s12, [r14, #68]    @ y[17] = x[12] from s12
vstr.32  s26, [r14, #72]    @ y[18] = x[26] from s26
vstr.32  s15, [r14, #76]    @ y[19] = x[15] from s15
str  r3, [r14, #80]    @ y[20] = x[2] from r3
vstr.32  s23, [r14, #84]    @ y[21] = x[23] from s23
str  r8, [r14, #88]    @ y[22] = x[6] from r8
vstr.32  s8, [r14, #92]    @ y[23] = x[8] from s8
vstr.32  s3, [r14, #96]    @ y[24] = x[3] from s3
str  r2, [r14, #100]    @ y[25] = x[20] from r2
str  r6, [r14, #104]    @ y[26] = x[30] from r6
vstr.32  s18, [r14, #108]    @ y[27] = x[18] from s18
vstr.32  s16, [r14, #112]    @ y[28] = x[16] from s16
str  r10, [r14, #116]    @ y[29] = x[10] from r10
vstr.32  s29, [r14, #120]    @ y[30] = x[29] from s29
vstr.32  s21, [r14, #124]    @ y[31] = x[21] from s21
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   decode_from_gft_full_length, .-decode_from_gft_full_length