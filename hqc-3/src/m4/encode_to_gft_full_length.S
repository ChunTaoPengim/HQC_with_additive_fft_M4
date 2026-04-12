.syntax unified
.cpu cortex-m4
.fpu fpv4-sp-d16
.global encode_to_gft_full_length
.type encode_to_gft_full_length, %function
.align 2
encode_to_gft_full_length:
push { r0-r12, r14 }
vpush {d8-d15}
mov  r14, r0
mov  r12, r1
vldm r12, {s0-s31}
vmov r0, s28    @ load x[28] from s28
vmov r1, s19    @ load x[19] from s19
eor r0, r0, r1    @ x[28] = x[28] ^ x[19]
vmov r2, s7    @ load x[7] from s7
vmov r3, s27    @ load x[27] from s27
eor r2, r2, r3    @ x[7] = x[7] ^ x[27]
vmov r4, s29    @ load x[29] from s29
vmov r5, s6    @ load x[6] from s6
eor r4, r4, r5    @ x[29] = x[29] ^ x[6]
vmov r6, s13    @ load x[13] from s13
vmov r7, s30    @ load x[30] from s30
eor r6, r6, r7    @ x[13] = x[13] ^ x[30]
vmov r8, s18    @ load x[18] from s18
eor r6, r6, r8    @ x[13] = x[13] ^ x[18]
vmov r9, s31    @ load x[31] from s31
vmov r10, s10    @ load x[10] from s10
eor r9, r9, r10    @ x[31] = x[31] ^ x[10]
vmov r11, s11    @ load x[11] from s11
vmov s10, r10    @ spill x[10] from r10
vmov r10, s24    @ load x[24] from s24
eor r11, r11, r10    @ x[11] = x[11] ^ x[24]
vmov s24, r10    @ spill x[24] from r10
vmov r10, s20    @ load x[20] from s20
eor r10, r10, r2    @ x[20] = x[20] ^ x[7]
vmov s7, r2    @ spill x[7] from r2
vmov r2, s26    @ load x[26] from s26
eor r7, r7, r2    @ x[30] = x[30] ^ x[26]
eor r10, r10, r6    @ x[20] = x[20] ^ x[13]
vmov s20, r10    @ spill x[20] from r10
vmov r10, s23    @ load x[23] from s23
eor r1, r1, r10    @ x[19] = x[19] ^ x[23]
vmov s23, r10    @ spill x[23] from r10
vmov r10, s1    @ load x[1] from s1
eor r1, r1, r10    @ x[19] = x[19] ^ x[1]
vmov s1, r10    @ spill x[1] from r10
vmov r10, s16    @ load x[16] from s16
vmov s13, r6    @ spill x[13] from r6
vmov r6, s12    @ load x[12] from s12
eor r10, r10, r6    @ x[16] = x[16] ^ x[12]
eor r8, r8, r6    @ x[18] = x[18] ^ x[12]
eor r6, r6, r4    @ x[12] = x[12] ^ x[29]
vmov s16, r10    @ spill x[16] from r10
vmov r10, s25    @ load x[25] from s25
vmov s18, r8    @ spill x[18] from r8
vmov r8, s5    @ load x[5] from s5
eor r10, r10, r8    @ x[25] = x[25] ^ x[5]
vmov s25, r10    @ spill x[25] from r10
vmov r10, s21    @ load x[21] from s21
eor r10, r10, r0    @ x[21] = x[21] ^ x[28]
eor r2, r2, r6    @ x[26] = x[26] ^ x[12]
vmov s12, r6    @ spill x[12] from r6
vmov r6, s4    @ load x[4] from s4
eor r8, r8, r6    @ x[5] = x[5] ^ x[4]
eor r1, r1, r7    @ x[19] = x[19] ^ x[30]
eor r4, r4, r6    @ x[29] = x[29] ^ x[4]
eor r7, r7, r4    @ x[30] = x[30] ^ x[29]
vmov s30, r7    @ spill x[30] from r7
vmov r7, s2    @ load x[2] from s2
vmov s4, r6    @ spill x[4] from r6
vmov r6, s8    @ load x[8] from s8
eor r7, r7, r6    @ x[2] = x[2] ^ x[8]
vmov s2, r7    @ spill x[2] from r7
vmov r7, s22    @ load x[22] from s22
eor r6, r6, r7    @ x[8] = x[8] ^ x[22]
vmov s22, r7    @ spill x[22] from r7
vmov r7, s15    @ load x[15] from s15
eor r7, r7, r9    @ x[15] = x[15] ^ x[31]
eor r8, r8, r3    @ x[5] = x[5] ^ x[27]
vmov s31, r9    @ spill x[31] from r9
vmov r9, s17    @ load x[17] from s17
eor r3, r3, r9    @ x[27] = x[27] ^ x[17]
vmov s17, r9    @ spill x[17] from r9
vmov r9, s9    @ load x[9] from s9
eor r8, r8, r9    @ x[5] = x[5] ^ x[9]
vmov s19, r1    @ spill x[19] from r1
vmov r1, s3    @ load x[3] from s3
eor r1, r1, r5    @ x[3] = x[3] ^ x[6]
eor r2, r2, r11    @ x[26] = x[26] ^ x[11]
vmov s11, r11    @ spill x[11] from r11
vmov r11, s18    @ load x[18] from s18
eor r11, r11, r10    @ x[18] = x[18] ^ x[21]
eor r11, r11, r6    @ x[18] = x[18] ^ x[8]
vmov s8, r6    @ spill x[8] from r6
vmov r6, s7    @ load x[7] from s7
eor r11, r11, r6    @ x[18] = x[18] ^ x[7]
vmov s18, r11    @ spill x[18] from r11
vmov r11, s13    @ load x[13] from s13
eor r11, r11, r1    @ x[13] = x[13] ^ x[3]
vmov s7, r6    @ spill x[7] from r6
vmov r6, s0    @ load x[0] from s0
eor r6, r6, r0    @ x[0] = x[0] ^ x[28]
vmov s28, r0    @ spill x[28] from r0
vmov r0, s12    @ load x[12] from s12
eor r0, r0, r5    @ x[12] = x[12] ^ x[6]
vmov s0, r6    @ spill x[0] from r6
vmov r6, s10    @ load x[10] from s10
eor r5, r5, r6    @ x[6] = x[6] ^ x[10]
eor r5, r5, r11    @ x[6] = x[6] ^ x[13]
eor r11, r11, r3    @ x[13] = x[13] ^ x[27]
eor r11, r11, r9    @ x[13] = x[13] ^ x[9]
eor r9, r9, r4    @ x[9] = x[9] ^ x[29]
eor r4, r4, r3    @ x[29] = x[29] ^ x[27]
eor r9, r9, r7    @ x[9] = x[9] ^ x[15]
eor r7, r7, r0    @ x[15] = x[15] ^ x[12]
vmov s13, r11    @ spill x[13] from r11
vmov r11, s24    @ load x[24] from s24
eor r0, r0, r11    @ x[12] = x[12] ^ x[24]
eor r3, r3, r11    @ x[27] = x[27] ^ x[24]
eor r7, r7, r10    @ x[15] = x[15] ^ x[21]
eor r10, r10, r11    @ x[21] = x[21] ^ x[24]
eor r11, r11, r6    @ x[24] = x[24] ^ x[10]
eor r6, r6, r8    @ x[10] = x[10] ^ x[5]
vmov s12, r0    @ spill x[12] from r0
vmov r0, s19    @ load x[19] from s19
eor r11, r11, r0    @ x[24] = x[24] ^ x[19]
eor r8, r8, r2    @ x[5] = x[5] ^ x[26]
vmov s29, r4    @ spill x[29] from r4
vmov r4, s31    @ load x[31] from s31
eor r2, r2, r4    @ x[26] = x[26] ^ x[31]
eor r0, r0, r4    @ x[19] = x[19] ^ x[31]
vmov s5, r8    @ spill x[5] from r8
vmov r8, s14    @ load x[14] from s14
eor r4, r4, r8    @ x[31] = x[31] ^ x[14]
eor r6, r6, r8    @ x[10] = x[10] ^ x[14]
vmov s10, r6    @ spill x[10] from r6
vmov r6, s16    @ load x[16] from s16
vmov s31, r4    @ spill x[31] from r4
vmov r4, s11    @ load x[11] from s11
eor r6, r6, r4    @ x[16] = x[16] ^ x[11]
vmov s19, r0    @ spill x[19] from r0
vmov r0, s20    @ load x[20] from s20
eor r8, r8, r0    @ x[14] = x[14] ^ x[20]
eor r2, r2, r0    @ x[26] = x[26] ^ x[20]
vmov s26, r2    @ spill x[26] from r2
vmov r2, s2    @ load x[2] from s2
eor r9, r9, r2    @ x[9] = x[9] ^ x[2]
eor r10, r10, r0    @ x[21] = x[21] ^ x[20]
vmov s2, r2    @ spill x[2] from r2
vmov r2, s22    @ load x[22] from s22
eor r10, r10, r2    @ x[21] = x[21] ^ x[22]
vmov s9, r9    @ spill x[9] from r9
vmov r9, s4    @ load x[4] from s4
eor r8, r8, r9    @ x[14] = x[14] ^ x[4]
vmov s14, r8    @ spill x[14] from r8
vmov r8, s1    @ load x[1] from s1
eor r0, r0, r8    @ x[20] = x[20] ^ x[1]
vmov s20, r0    @ spill x[20] from r0
vmov r0, s30    @ load x[30] from s30
eor r2, r2, r0    @ x[22] = x[22] ^ x[30]
eor r8, r8, r1    @ x[1] = x[1] ^ x[3]
vmov s22, r2    @ spill x[22] from r2
vmov r2, s17    @ load x[17] from s17
eor r11, r11, r2    @ x[24] = x[24] ^ x[17]
vmov s24, r11    @ spill x[24] from r11
vmov r11, s7    @ load x[7] from s7
eor r0, r0, r11    @ x[30] = x[30] ^ x[7]
eor r11, r11, r4    @ x[7] = x[7] ^ x[11]
eor r10, r10, r4    @ x[21] = x[21] ^ x[11]
eor r1, r1, r3    @ x[3] = x[3] ^ x[27]
vmov s21, r10    @ spill x[21] from r10
vmov r10, s23    @ load x[23] from s23
eor r3, r3, r10    @ x[27] = x[27] ^ x[23]
eor r7, r7, r5    @ x[15] = x[15] ^ x[6]
eor r2, r2, r10    @ x[17] = x[17] ^ x[23]
vmov s27, r3    @ spill x[27] from r3
vmov r3, s8    @ load x[8] from s8
eor r3, r3, r10    @ x[8] = x[8] ^ x[23]
eor r7, r7, r8    @ x[15] = x[15] ^ x[1]
eor r8, r8, r9    @ x[1] = x[1] ^ x[4]
eor r9, r9, r0    @ x[4] = x[4] ^ x[30]
eor r10, r10, r11    @ x[23] = x[23] ^ x[7]
vmov s1, r8    @ spill x[1] from r8
vmov r8, s19    @ load x[19] from s19
eor r6, r6, r8    @ x[16] = x[16] ^ x[19]
vmov s15, r7    @ spill x[15] from r7
vmov r7, s14    @ load x[14] from s14
eor r0, r0, r7    @ x[30] = x[30] ^ x[14]
vmov s30, r0    @ spill x[30] from r0
vmov r0, s24    @ load x[24] from s24
vmov s16, r6    @ spill x[16] from r6
vmov r6, s5    @ load x[5] from s5
eor r0, r0, r6    @ x[24] = x[24] ^ x[5]
vmov s24, r0    @ spill x[24] from r0
vmov r0, s0    @ load x[0] from s0
vmov s5, r6    @ spill x[5] from r6
vmov r6, s9    @ load x[9] from s9
eor r0, r0, r6    @ x[0] = x[0] ^ x[9]
vmov s9, r6    @ spill x[9] from r6
vmov r6, s18    @ load x[18] from s18
eor r0, r0, r6    @ x[0] = x[0] ^ x[18]
vmov s0, r0    @ spill x[0] from r0
vmov r0, s20    @ load x[20] from s20
eor r4, r4, r0    @ x[11] = x[11] ^ x[20]
eor r0, r0, r2    @ x[20] = x[20] ^ x[17]
vmov s14, r7    @ spill x[14] from r7
vmov r7, s29    @ load x[29] from s29
eor r10, r10, r7    @ x[23] = x[23] ^ x[29]
eor r11, r11, r3    @ x[7] = x[7] ^ x[8]
vmov s23, r10    @ spill x[23] from r10
vmov r10, s2    @ load x[2] from s2
eor r3, r3, r10    @ x[8] = x[8] ^ x[2]
eor r7, r7, r0    @ x[29] = x[29] ^ x[20]
eor r0, r0, r9    @ x[20] = x[20] ^ x[4]
eor r9, r9, r10    @ x[4] = x[4] ^ x[2]
eor r10, r10, r2    @ x[2] = x[2] ^ x[17]
vmov s8, r3    @ spill x[8] from r3
vmov r3, s25    @ load x[25] from s25
eor r2, r2, r3    @ x[17] = x[17] ^ x[25]
eor r9, r9, r8    @ x[4] = x[4] ^ x[19]
vmov s19, r8    @ spill x[19] from r8
vmov r8, s31    @ load x[31] from s31
eor r0, r0, r8    @ x[20] = x[20] ^ x[31]
eor r8, r8, r6    @ x[31] = x[31] ^ x[18]
eor r10, r10, r6    @ x[2] = x[2] ^ x[18]
vmov s20, r0    @ spill x[20] from r0
vmov r0, s26    @ load x[26] from s26
eor r4, r4, r0    @ x[11] = x[11] ^ x[26]
eor r11, r11, r1    @ x[7] = x[7] ^ x[3]
vmov s2, r10    @ spill x[2] from r10
vmov r10, s28    @ load x[28] from s28
eor r11, r11, r10    @ x[7] = x[7] ^ x[28]
eor r4, r4, r9    @ x[11] = x[11] ^ x[4]
vmov s7, r11    @ spill x[7] from r11
vmov r11, s12    @ load x[12] from s12
eor r6, r6, r11    @ x[18] = x[18] ^ x[12]
eor r11, r11, r2    @ x[12] = x[12] ^ x[17]
vmov s11, r4    @ spill x[11] from r4
vmov r4, s10    @ load x[10] from s10
eor r7, r7, r4    @ x[29] = x[29] ^ x[10]
eor r11, r11, r5    @ x[12] = x[12] ^ x[6]
vmov s29, r7    @ spill x[29] from r7
vmov r7, s19    @ load x[19] from s19
eor r7, r7, r3    @ x[19] = x[19] ^ x[25]
eor r3, r3, r5    @ x[25] = x[25] ^ x[6]
eor r9, r9, r5    @ x[4] = x[4] ^ x[6]
eor r3, r3, r8    @ x[25] = x[25] ^ x[31]
vmov s25, r3    @ spill x[25] from r3
vmov r3, s14    @ load x[14] from s14
eor r8, r8, r3    @ x[31] = x[31] ^ x[14]
eor r3, r3, r2    @ x[14] = x[14] ^ x[17]
eor r2, r2, r10    @ x[17] = x[17] ^ x[28]
vmov s14, r3    @ spill x[14] from r3
vmov r3, s5    @ load x[5] from s5
eor r5, r5, r3    @ x[6] = x[6] ^ x[5]
vmov s4, r9    @ spill x[4] from r9
vmov r9, s21    @ load x[21] from s21
vmov s10, r4    @ spill x[10] from r4
vmov r4, s27    @ load x[27] from s27
eor r9, r9, r4    @ x[21] = x[21] ^ x[27]
eor r10, r10, r3    @ x[28] = x[28] ^ x[5]
vmov s21, r9    @ spill x[21] from r9
vmov r9, s13    @ load x[13] from s13
eor r2, r2, r9    @ x[17] = x[17] ^ x[13]
vmov s13, r9    @ spill x[13] from r9
vmov r9, s8    @ load x[8] from s8
eor r3, r3, r9    @ x[5] = x[5] ^ x[8]
vmov s5, r3    @ spill x[5] from r3
vmov r3, s16    @ load x[16] from s16
eor r4, r4, r3    @ x[27] = x[27] ^ x[16]
eor r11, r11, r0    @ x[12] = x[12] ^ x[26]
vmov s12, r11    @ spill x[12] from r11
vmov r11, s22    @ load x[22] from s22
eor r5, r5, r11    @ x[6] = x[6] ^ x[22]
eor r4, r4, r11    @ x[27] = x[27] ^ x[22]
eor r8, r8, r10    @ x[31] = x[31] ^ x[28]
vmov s27, r4    @ spill x[27] from r4
vmov r4, s24    @ load x[24] from s24
eor r10, r10, r4    @ x[28] = x[28] ^ x[24]
eor r8, r8, r9    @ x[31] = x[31] ^ x[8]
vmov s28, r10    @ spill x[28] from r10
vmov r10, s9    @ load x[9] from s9
eor r0, r0, r10    @ x[26] = x[26] ^ x[9]
vmov s26, r0    @ spill x[26] from r0
vmov r0, s15    @ load x[15] from s15
eor r4, r4, r0    @ x[24] = x[24] ^ x[15]
eor r1, r1, r4    @ x[3] = x[3] ^ x[24]
eor r8, r8, r1    @ x[31] = x[31] ^ x[3]
eor r7, r7, r0    @ x[19] = x[19] ^ x[15]
vmov s31, r8    @ spill x[31] from r8
vmov r8, s23    @ load x[23] from s23
eor r6, r6, r8    @ x[18] = x[18] ^ x[23]
eor r5, r5, r2    @ x[6] = x[6] ^ x[17]
eor r4, r4, r8    @ x[24] = x[24] ^ x[23]
eor r4, r4, r11    @ x[24] = x[24] ^ x[22]
vmov s6, r5    @ spill x[6] from r5
vmov r5, s1    @ load x[1] from s1
eor r8, r8, r5    @ x[23] = x[23] ^ x[1]
eor r11, r11, r5    @ x[22] = x[22] ^ x[1]
eor r5, r5, r1    @ x[1] = x[1] ^ x[3]
eor r9, r9, r6    @ x[8] = x[8] ^ x[18]
eor r1, r1, r6    @ x[3] = x[3] ^ x[18]
vmov s22, r11    @ spill x[22] from r11
vmov r11, s10    @ load x[10] from s10
eor r6, r6, r11    @ x[18] = x[18] ^ x[10]
vmov s19, r7    @ spill x[19] from r7
vmov r7, s11    @ load x[11] from s11
eor r7, r7, r9    @ x[11] = x[11] ^ x[8]
vmov s11, r7    @ spill x[11] from r7
vmov r7, s13    @ load x[13] from s13
eor r11, r11, r7    @ x[10] = x[10] ^ x[13]
eor r9, r9, r10    @ x[8] = x[8] ^ x[9]
vmov s8, r9    @ spill x[8] from r9
vmov r9, s4    @ load x[4] from s4
eor r8, r8, r9    @ x[23] = x[23] ^ x[4]
vmov s3, r1    @ spill x[3] from r1
vmov r1, s21    @ load x[21] from s21
eor r2, r2, r1    @ x[17] = x[17] ^ x[21]
eor r11, r11, r10    @ x[10] = x[10] ^ x[9]
eor r5, r5, r3    @ x[1] = x[1] ^ x[16]
vmov s21, r1    @ spill x[21] from r1
vmov r1, s7    @ load x[7] from s7
eor r6, r6, r1    @ x[18] = x[18] ^ x[7]
eor r7, r7, r9    @ x[13] = x[13] ^ x[4]
eor r1, r1, r8    @ x[7] = x[7] ^ x[23]
eor r9, r9, r1    @ x[4] = x[4] ^ x[7]
eor r1, r1, r3    @ x[7] = x[7] ^ x[16]
vmov s18, r6    @ spill x[18] from r6
vmov r6, s30    @ load x[30] from s30
eor r6, r6, r0    @ x[30] = x[30] ^ x[15]
vmov s7, r1    @ spill x[7] from r1
vmov r1, s3    @ load x[3] from s3
eor r1, r1, r4    @ x[3] = x[3] ^ x[24]
vmov s30, r6    @ spill x[30] from r6
vmov r6, s29    @ load x[29] from s29
eor r0, r0, r6    @ x[15] = x[15] ^ x[29]
vmov s15, r0    @ spill x[15] from r0
vmov r0, s19    @ load x[19] from s19
eor r3, r3, r0    @ x[16] = x[16] ^ x[19]
vmov s3, r1    @ spill x[3] from r1
vmov r1, s2    @ load x[2] from s2
vmov s19, r0    @ spill x[19] from r0
vmov r0, s5    @ load x[5] from s5
eor r1, r1, r0    @ x[2] = x[2] ^ x[5]
vmov s5, r0    @ spill x[5] from r0
vmov r0, s28    @ load x[28] from s28
eor r9, r9, r0    @ x[4] = x[4] ^ x[28]
vmov s29, r6    @ spill x[29] from r6
vmov r6, s22    @ load x[22] from s22
eor r1, r1, r6    @ x[2] = x[2] ^ x[22]
vmov s22, r6    @ spill x[22] from r6
vmov r6, s0    @ load x[0] from s0
vmov s23, r8    @ spill x[23] from r8
vmov r8, s25    @ load x[25] from s25
eor r6, r6, r8    @ x[0] = x[0] ^ x[25]
vmov s10, r11    @ spill x[10] from r11
vmov r11, s27    @ load x[27] from s27
eor r7, r7, r11    @ x[13] = x[13] ^ x[27]
vmov s13, r7    @ spill x[13] from r7
vmov r7, s26    @ load x[26] from s26
eor r11, r11, r7    @ x[27] = x[27] ^ x[26]
vmov s27, r11    @ spill x[27] from r11
vmov r11, s20    @ load x[20] from s20
eor r11, r11, r0    @ x[20] = x[20] ^ x[28]
eor r10, r10, r3    @ x[9] = x[9] ^ x[16]
eor r4, r4, r1    @ x[24] = x[24] ^ x[2]
vmov s16, r3    @ spill x[16] from r3
vmov r3, s11    @ load x[11] from s11
eor r3, r3, r5    @ x[11] = x[11] ^ x[1]
eor r5, r5, r8    @ x[1] = x[1] ^ x[25]
eor r4, r4, r2    @ x[24] = x[24] ^ x[17]
eor r6, r6, r9    @ x[0] = x[0] ^ x[4]
vmov s24, r4    @ spill x[24] from r4
vmov r4, s10    @ load x[10] from s10
eor r9, r9, r4    @ x[4] = x[4] ^ x[10]
vmov s9, r10    @ spill x[9] from r10
vmov r10, s31    @ load x[31] from s31
eor r10, r10, r7    @ x[31] = x[31] ^ x[26]
vmov s1, r5    @ spill x[1] from r5
vmov r5, s23    @ load x[23] from s23
eor r4, r4, r5    @ x[10] = x[10] ^ x[23]
eor r5, r5, r2    @ x[23] = x[23] ^ x[17]
vmov s17, r2    @ spill x[17] from r2
vmov r2, s29    @ load x[29] from s29
vmov s4, r9    @ spill x[4] from r9
vmov r9, s6    @ load x[6] from s6
eor r2, r2, r9    @ x[29] = x[29] ^ x[6]
eor r0, r0, r2    @ x[28] = x[28] ^ x[29]
str  r6, [r14, #0]    @ y[0] = x[0] from r6
str  r7, [r14, #4]    @ y[1] = x[26] from r7
vstr.32  s7, [r14, #8]    @ y[2] = x[7] from s7
str  r5, [r14, #12]    @ y[3] = x[23] from r5
vstr.32  s24, [r14, #16]    @ y[4] = x[24] from s24
vstr.32  s27, [r14, #20]    @ y[5] = x[27] from s27
vstr.32  s22, [r14, #24]    @ y[6] = x[22] from s22
vstr.32  s13, [r14, #28]    @ y[7] = x[13] from s13
str  r4, [r14, #32]    @ y[8] = x[10] from r4
vstr.32  s9, [r14, #36]    @ y[9] = x[9] from s9
str  r1, [r14, #40]    @ y[10] = x[2] from r1
vstr.32  s1, [r14, #44]    @ y[11] = x[1] from s1
vstr.32  s16, [r14, #48]    @ y[12] = x[16] from s16
vstr.32  s30, [r14, #52]    @ y[13] = x[30] from s30
vstr.32  s5, [r14, #56]    @ y[14] = x[5] from s5
vstr.32  s3, [r14, #60]    @ y[15] = x[3] from s3
str  r3, [r14, #64]    @ y[16] = x[11] from r3
str  r8, [r14, #68]    @ y[17] = x[25] from r8
vstr.32  s14, [r14, #72]    @ y[18] = x[14] from s14
vstr.32  s4, [r14, #76]    @ y[19] = x[4] from s4
str  r2, [r14, #80]    @ y[20] = x[29] from r2
vstr.32  s18, [r14, #84]    @ y[21] = x[18] from s18
str  r9, [r14, #88]    @ y[22] = x[6] from r9
vstr.32  s8, [r14, #92]    @ y[23] = x[8] from s8
str  r11, [r14, #96]    @ y[24] = x[20] from r11
str  r0, [r14, #100]    @ y[25] = x[28] from r0
vstr.32  s17, [r14, #104]    @ y[26] = x[17] from s17
vstr.32  s12, [r14, #108]    @ y[27] = x[12] from s12
vstr.32  s15, [r14, #112]    @ y[28] = x[15] from s15
str  r10, [r14, #116]    @ y[29] = x[31] from r10
vstr.32  s19, [r14, #120]    @ y[30] = x[19] from s19
vstr.32  s21, [r14, #124]    @ y[31] = x[21] from s21
vpop { d8-d15 }
pop { r0-r12, r14 }
bx lr
.size   encode_to_gft_full_length, .-encode_to_gft_full_length