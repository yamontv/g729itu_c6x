/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder with Annex B    ANSI-C Source Code
 Version 1.3    Last modified: August 1997

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies,
 Rockwell International
 All rights reserved.
 */

extern const Word16 hamwindow[L_WINDOW];
extern const Word16 lag_h[M + 2];
extern const Word16 lag_l[M + 2];
extern const Word16 table[65];
extern const Word16 slope[64];
extern const Word16 table2[64];
extern const Word16 slope_cos[64];
extern const Word16 slope_acos[64];
extern const Word16 lspcb1[NC0][M];
extern const Word16 lspcb2[NC1][M];
extern const Word16 fg[2][MA_NP][M];
extern const Word16 fg_sum[2][M];
extern const Word16 fg_sum_inv[2][M];
extern const Word16 grid[GRID_POINTS + 1];
extern const Word16 freq_prev_reset[M];
extern const Word16 inter_3l[FIR_SIZE_SYN];
extern const Word16 pred[4];
extern const Word16 gbk1[NCODE1][2];
extern const Word16 gbk2[NCODE2][2];
extern const Word16 map1[NCODE1];
extern const Word16 map2[NCODE2];
extern const Word16 coef[2][2];
extern const Word32 L_coef[2][2];
extern const Word16 thr1[NCODE1 - NCAN1];
extern const Word16 thr2[NCODE2 - NCAN2];
extern const Word16 imap1[NCODE1];
extern const Word16 imap2[NCODE2];
extern const Word16 b100[3];
extern const Word16 a100[3];
extern const Word16 b140[3];
extern const Word16 a140[3];
extern const Word16 bitsno[PRM_SIZE];
extern const Word16 bitsno2[4];
extern const Word16 tabpow[33];
extern const Word16 tablog[33];
extern const Word16 tabsqr[49];
extern const Word16 tab_zone[PIT_MAX + L_INTERPOL - 1];

