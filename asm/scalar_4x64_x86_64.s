/***********************************************************************
 * Copyright (c) 2013, 2014 Pieter Wuille                              *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php *
 ***********************************************************************/
/*
 * x86_64 scalar_4x64 assembly.
 * Extracted from libsecp256k1 src/scalar_4x64_impl.h (USE_ASM_X86_64).
 * Assembly instructions are copied verbatim from libsecp256k1.
 */

	.text

/* SECP256K1_N_C_0 = ~N_0 + 1 = 0x402DA1732FC9BEBF */
/* SECP256K1_N_C_1 = ~N_1     = 0x4551231950B75FC4 */

/*
 * blvm_secp256k1_scalar_mul_512(l8, a, b)
 * SysV: rdi=l8, rsi=a, rdx=b
 * libsecp256k1 asm expects: rdi=a, rsi=l8, rdx=b
 * Swap rdi and rsi at entry.
 */
	.globl	blvm_secp256k1_scalar_mul_512
	.hidden	blvm_secp256k1_scalar_mul_512
	.type	blvm_secp256k1_scalar_mul_512, @function
blvm_secp256k1_scalar_mul_512:
	pushq	%rbx
	pushq	%r12
	pushq	%r13
	pushq	%r14
	pushq	%r15
	xchgq	%rdi, %rsi
	/* Preload - from libsecp256k1 scalar_4x64_impl.h lines 682-688 */
	movq	0(%rdi), %r15
	movq	8(%rdi), %rbx
	movq	16(%rdi), %rcx
	movq	0(%rdx), %r11
	movq	8(%rdx), %r12
	movq	16(%rdx), %r13
	movq	24(%rdx), %r14
	/* (rax,rdx) = a0 * b0 */
	movq	%r15, %rax
	mulq	%r11
	/* Extract l8[0] */
	movq	%rax, 0(%rsi)
	/* (r8,r9,r10) = (rdx) */
	movq	%rdx, %r8
	xorq	%r9, %r9
	xorq	%r10, %r10
	/* (r8,r9,r10) += a0 * b1 */
	movq	%r15, %rax
	mulq	%r12
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	/* (r8,r9,r10) += a1 * b0 */
	movq	%rbx, %rax
	mulq	%r11
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	/* Extract l8[1] */
	movq	%r8, 8(%rsi)
	xorq	%r8, %r8
	/* (r9,r10,r8) += a0 * b2 */
	movq	%r15, %rax
	mulq	%r13
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	/* (r9,r10,r8) += a1 * b1 */
	movq	%rbx, %rax
	mulq	%r12
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	/* (r9,r10,r8) += a2 * b0 */
	movq	%rcx, %rax
	mulq	%r11
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	/* Extract l8[2] */
	movq	%r9, 16(%rsi)
	xorq	%r9, %r9
	/* (r10,r8,r9) += a0 * b3 */
	movq	%r15, %rax
	mulq	%r14
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	/* Preload a3 */
	movq	24(%rdi), %r15
	/* (r10,r8,r9) += a1 * b2 */
	movq	%rbx, %rax
	mulq	%r13
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	/* (r10,r8,r9) += a2 * b1 */
	movq	%rcx, %rax
	mulq	%r12
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	/* (r10,r8,r9) += a3 * b0 */
	movq	%r15, %rax
	mulq	%r11
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	/* Extract l8[3] */
	movq	%r10, 24(%rsi)
	xorq	%r10, %r10
	/* (r8,r9,r10) += a1 * b3 */
	movq	%rbx, %rax
	mulq	%r14
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	/* (r8,r9,r10) += a2 * b2 */
	movq	%rcx, %rax
	mulq	%r13
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	/* (r8,r9,r10) += a3 * b1 */
	movq	%r15, %rax
	mulq	%r12
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	/* Extract l8[4] */
	movq	%r8, 32(%rsi)
	xorq	%r8, %r8
	/* (r9,r10,r8) += a2 * b3 */
	movq	%rcx, %rax
	mulq	%r14
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	/* (r9,r10,r8) += a3 * b2 */
	movq	%r15, %rax
	mulq	%r13
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	/* Extract l8[5] */
	movq	%r9, 40(%rsi)
	/* (r10,r8) += a3 * b3 */
	movq	%r15, %rax
	mulq	%r14
	addq	%rax, %r10
	adcq	%rdx, %r8
	/* Extract l8[6] */
	movq	%r10, 48(%rsi)
	/* Extract l8[7] */
	movq	%r8, 56(%rsi)
	popq	%r15
	popq	%r14
	popq	%r13
	popq	%r12
	popq	%rbx
	ret

/*
 * blvm_secp256k1_scalar_reduce_512(r, l) -> returns overflow (rax)
 * SysV: rdi=r, rsi=l
 * Assembly from libsecp256k1 scalar_4x64_impl.h (three blocks).
 */
	.globl	blvm_secp256k1_scalar_reduce_512
	.hidden	blvm_secp256k1_scalar_reduce_512
	.type	blvm_secp256k1_scalar_reduce_512, @function
blvm_secp256k1_scalar_reduce_512:
	pushq	%r12
	pushq	%r13
	pushq	%r14
	subq	$96, %rsp
	/* Block 1: Reduce 512 bits into 385. rsi=l. Output m0..m6 to stack. */
	movq	32(%rsi), %r11
	movq	40(%rsi), %r12
	movq	48(%rsi), %r13
	movq	56(%rsi), %r14
	movq	0(%rsi), %r8
	xorq	%r9, %r9
	xorq	%r10, %r10
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r11
	addq	%rax, %r8
	adcq	%rdx, %r9
	movq	%r8, 0(%rsp)
	xorq	%r8, %r8
	addq	8(%rsi), %r9
	adcq	$0, %r10
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r12
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r11
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	movq	%r9, 8(%rsp)
	xorq	%r9, %r9
	addq	16(%rsi), %r10
	adcq	$0, %r8
	adcq	$0, %r9
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r13
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r12
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	addq	%r11, %r10
	adcq	$0, %r8
	adcq	$0, %r9
	movq	%r10, 16(%rsp)
	xorq	%r10, %r10
	addq	24(%rsi), %r8
	adcq	$0, %r9
	adcq	$0, %r10
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r14
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r13
	addq	%rax, %r8
	adcq	%rdx, %r9
	adcq	$0, %r10
	addq	%r12, %r8
	adcq	$0, %r9
	adcq	$0, %r10
	movq	%r8, 24(%rsp)
	xorq	%r8, %r8
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r14
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	addq	%r13, %r9
	adcq	$0, %r10
	adcq	$0, %r8
	movq	%r9, 32(%rsp)
	addq	%r14, %r10
	adcq	$0, %r8
	movq	%r10, 40(%rsp)
	movq	%r8, 48(%rsp)

	/* Block 2: Reduce 385 bits into 258. m0..m6 from stack, output p0..p4 to stack. */
	movq	32(%rsp), %r11
	movq	40(%rsp), %r12
	movq	48(%rsp), %r13
	movq	0(%rsp), %r8
	xorq	%r9, %r9
	xorq	%r10, %r10
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r11
	addq	%rax, %r8
	adcq	%rdx, %r9
	movq	%r8, 56(%rsp)
	xorq	%r8, %r8
	addq	8(%rsp), %r9
	adcq	$0, %r10
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r12
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r11
	addq	%rax, %r9
	adcq	%rdx, %r10
	adcq	$0, %r8
	movq	%r9, 64(%rsp)
	xorq	%r9, %r9
	addq	16(%rsp), %r10
	adcq	$0, %r8
	adcq	$0, %r9
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r13
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r12
	addq	%rax, %r10
	adcq	%rdx, %r8
	adcq	$0, %r9
	addq	%r11, %r10
	adcq	$0, %r8
	adcq	$0, %r9
	movq	%r10, 72(%rsp)
	addq	24(%rsp), %r8
	adcq	$0, %r9
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r13
	addq	%rax, %r8
	adcq	%rdx, %r9
	addq	%r12, %r8
	adcq	$0, %r9
	movq	%r8, 80(%rsp)
	addq	%r13, %r9
	movq	%r9, 88(%rsp)

	/* Block 3: Reduce 258 bits into 256. p0..p4 from stack, r in rdi. Return c in rax. */
	movq	88(%rsp), %r10
	movabsq	$0x402DA1732FC9BEBF, %rax
	mulq	%r10
	addq	56(%rsp), %rax
	adcq	$0, %rdx
	movq	%rax, 0(%rdi)
	movq	%rdx, %r8
	xorq	%r9, %r9
	addq	64(%rsp), %r8
	adcq	$0, %r9
	movabsq	$0x4551231950B75FC4, %rax
	mulq	%r10
	addq	%rax, %r8
	adcq	%rdx, %r9
	movq	%r8, 8(%rdi)
	xorq	%r8, %r8
	addq	%r10, %r9
	adcq	$0, %r8
	addq	72(%rsp), %r9
	adcq	$0, %r8
	movq	%r9, 16(%rdi)
	xorq	%r9, %r9
	addq	80(%rsp), %r8
	adcq	$0, %r9
	movq	%r8, 24(%rdi)
	movq	%r9, %rax

	addq	$96, %rsp
	popq	%r14
	popq	%r13
	popq	%r12
	ret
