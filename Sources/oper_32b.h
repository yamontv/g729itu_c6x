/*
 ITU-T G.729A Speech Coder    ANSI-C Source Code
 Version 1.1    Last modified: September 1996

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies
 All rights reserved.
 */

/* Double precision operations */

static inline void L_Extract(int L_32, short *hi, short *lo)
{
    *hi = extract_h(L_32);
    *lo = (L_32 & 0xFFFF) >> 1;  // lo = L_32>>1
    return;
}

#define L_Comp(hi, lo)  (((Word32)(hi) << 16) + ((Word32)(lo) << 1))

static inline int Mpy_32(short hi1, short lo1, short hi2, short lo2)
{
    int L_32;

    L_32 = L_mult(hi1, hi2);
    L_32 = L_mac(L_32, mult(hi1, lo2), 1);
    L_32 = L_mac(L_32, mult(lo1, hi2), 1);

    return (L_32);
}

static inline Word32 Mpy_32_16(Word16 hi, Word16 lo, Word16 n)
{
    Word32 L_32;

    L_32 = L_mult(hi, n);
    L_32 = L_mac(L_32, mult(lo, n), 1);

    return (L_32);
}

static inline Word32 Div_32(Word32 L_num, Word16 denom_hi, Word16 denom_lo)
{
    Word16 approx, hi, lo, n_hi, n_lo;
    Word32 L_32;

    /* First approximation: 1 / L_denom = 1/denom_hi */

    approx = div_s((Word16) 0x3fff, denom_hi); /* result in Q14 */
    /* Note: 3fff = 0.5 in Q15 */

    /* 1/L_denom = approx * (2.0 - L_denom * approx) */

    L_32 = Mpy_32_16(denom_hi, denom_lo, approx); /* result in Q30 */

    L_32 = L_sub(MAX_32, L_32); /* result in Q30 */

    hi = (Word16)(L_32 >> 16);
    lo = (L_32 >> 1) - (hi << 15);

    L_32 = Mpy_32_16(hi, lo, approx); /* = 1/L_denom in Q29 */

    /* L_num * (1/L_denom) */
    hi = (Word16)(L_32 >> 16);
    lo = (L_32 >> 1) - (hi << 15);
    n_hi = (Word16)(L_num >> 16);
    n_lo = (L_num >> 1) - (n_hi << 15);
    L_32 = Mpy_32(n_hi, n_lo, hi, lo); /* result in Q29   */
    L_32 = L_shl(L_32, 2); /* From Q29 to Q31 */

    return (L_32);
}
