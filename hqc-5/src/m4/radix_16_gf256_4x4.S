.syntax unified
.cpu cortex-m4
.thumb
.fpu fpv4-sp-d16
.global radix_16_gf256_4x4
.type radix_16_gf256_4x4, %function
.align 2
radix_16_gf256_4x4:

stmdb   sp!, { r4-r11, lr }

@ mask32 = 0x11111111
mov.w	r2, #0x11111111
@ _0x1D_RDX16_ = 0x00011101
movw    r9, #0x1101
movt    r9, #0x0001

@ a0 = (a)&mask32;
and.w	r3, r2, r0
@ b0 = (b)&mask32;
and.w	r4, r2, r1
@ r6||r5 = a0*b0
umull   r5, r6, r3, r4
@ r6 &= mask32
and.w   r6, r6, r2
@ uint64_t r1 = rdx16_mul( x_hi8bit , _0x1D_RDX16_ );
umull   r12, r10, r6, r9
@ uint32_t r2l = mla( r1h, _0x1D_RDX16_ , r1l );  // multiplication without normalization
mla     r12, r10, r9, r12
@ (((uint32_t)t0) ^ rdx16_gf256_reduce( t0>>32 )) = ((uint32_t)t0) ^ r10
eor     r5, r5, r12
@ r2l &= _NORM_RDX16_;
and.w   r5, r5, r2

@ a2 = (a >> 1)&mask32;
and.w	r3, r2, r0, lsr#1
@ b2 = (b >> 1)&mask32;
and.w	r4, r2, r1, lsr#1
@ r7||r6 = a1*b1
umull   r6, r7, r3, r4
@ r7 &= mask32
and.w   r7, r7, r2
@ uint64_t r1 = rdx16_mul( x_hi8bit , _0x1D_RDX16_ );
umull   r12, r10, r7, r9
@ uint32_t r2l = mla( r1h, _0x1D_RDX16_ , r1l );  // multiplication without normalization
mla     r12, r10, r9, r12
@ (((uint32_t)t0) ^ rdx16_gf256_reduce( t0>>32 )) = ((uint32_t)t0) ^ r10
eor     r6, r6, r12
@ r2l &= _NORM_RDX16_;
and.w   r6, r6, r2

@ a1 = (a >> 2)&mask32;
and.w	r3, r2, r0, lsr#2
@ b1 = (b >> 2)&mask32;
and.w	r4, r2, r1, lsr#2
@ r8||r7 = a1*b1
umull   r7, r8, r3, r4
@ r8 &= mask32
and.w   r8, r8, r2
@ uint64_t r1 = rdx16_mul( x_hi8bit , _0x1D_RDX16_ );
umull   r12, r10, r8, r9
@ uint32_t r2l = mla( r1h, _0x1D_RDX16_ , r1l );  // multiplication without normalization
mla     r12, r10, r9, r12
@ (((uint32_t)t0) ^ rdx16_gf256_reduce( t0>>32 )) = ((uint32_t)t0) ^ r10
eor     r7, r7, r12
@ r2l &= _NORM_RDX16_;
and.w   r7, r7, r2

@ a3 = (a >> 3)&mask32;
and.w	r3, r2, r0, lsr#3
@ b3 = (b >> 3)&mask32;
and.w	r4, r2, r1, lsr#3
@ r1||r8 = a3*b3
umull   r8, r1, r3, r4
@ r1 &= mask32
and.w   r1, r1, r2
@ uint64_t r1 = rdx16_mul( x_hi8bit , _0x1D_RDX16_ );
umull   r12, r10, r1, r9
@ uint32_t r2l = mla( r1h, _0x1D_RDX16_ , r1l );  // multiplication without normalization
mla     r12, r10, r9, r12
@ (((uint32_t)t0) ^ rdx16_gf256_reduce( t0>>32 )) = ((uint32_t)t0) ^ r10
eor     r8, r8, r12
@ r2l &= _NORM_RDX16_;
and.w   r8, r8, r2

orr.w   r5, r5, r6, lsl#1
orr.w   r5, r5, r7, lsl#2
orr.w   r0, r5, r8, lsl#3

ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}
