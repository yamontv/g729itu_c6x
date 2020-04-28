/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Annex B     ANSI-C Source Code
 Version 1.3    Last modified: August 1997
 Copyright (c) 1996, France Telecom, Rockwell International,
 Universite de Sherbrooke.
 All rights reserved.
 */

/* DTX and Comfort Noise Generator - Encoder part */

#include <stdio.h>
#include <stdlib.h>

#include "typedef.h"
#include "basic_op.h"
#include "ld8a.h"
#include "oper_32b.h"
#include "tab_ld8a.h"
#include "tab_dtx.h"
#include "sid.h"

#define FR_SID_MIN      3
#define FRAC_THRESH1    4855
#define FRAC_THRESH2    3161

/* Local functions */
static void Calc_pastfilt(struct dtx_ctx *dtxc);
static void Calc_RCoeff(Word16 *Coeff, Word16 *RCoeff, Word16 *sh_RCoeff);
static Word16 Cmp_filt(Word16 *RCoeff, Word16 sh_RCoeff, Word16 *acf,
		Word16 alpha, Word16 Fracthresh);
static void Calc_sum_acf(Word16 *acf, Word16 *sh_acf, Word16 *sum,
		Word16 *sh_sum, Word16 nb);
static void Update_sumAcf(struct dtx_ctx *dtxc);

/* Last A(z) for case of unstable filter */
static const Word16 old_A_c[M + 1] = { 4096, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static const Word16 old_rc_c[2] = { 0, 0 };
/*-----------------------------------------------------------*
 * procedure Init_Cod_cng:                                   *
 *           ~~~~~~~~~~~~                                    *
 *   Initialize variables used for dtx at the encoder        *
 *-----------------------------------------------------------*/
void Init_Cod_cng(struct dtx_ctx *dtxc) {
	Word16 i;

	for (i = 0; i < SIZ_SUMACF; i++)
		dtxc->sumAcf[i] = 0;
	for (i = 0; i < NB_SUMACF; i++)
		dtxc->sh_sumAcf[i] = 40;

	for (i = 0; i < SIZ_ACF; i++)
		dtxc->Acf[i] = 0;
	for (i = 0; i < NB_CURACF; i++)
		dtxc->sh_Acf[i] = 40;

	for (i = 0; i < NB_GAIN; i++)
		dtxc->sh_ener[i] = 40;
	for (i = 0; i < NB_GAIN; i++)
		dtxc->ener[i] = 0;

	dtxc->cur_gain = 0;
	dtxc->fr_cur = 0;
	dtxc->flag_chang = 0;

	Copy(old_A_c, dtxc->old_A, M + 1);
	Copy(old_rc_c, dtxc->old_rc, 2);

	return;
}

/*-----------------------------------------------------------*
 * procedure Cod_cng:                                        *
 *           ~~~~~~~~                                        *
 *   computes DTX decision                                   *
 *   encodes SID frames                                      *
 *   computes CNG excitation for encoder update              *
 *-----------------------------------------------------------*/
void Cod_cng(struct dtx_ctx * restrict dtxc, Word16 *exc, /* (i/o) : excitation array                     */
Word16 pastVad, /* (i)   : previous VAD decision                */
Word16 * restrict lsp_old_q, /* (i/o) : previous quantized lsp               */
Word16 *Aq, /* (o)   : set of interpolated LPC coefficients */
Word16 *ana, /* (o)   : coded SID parameters                 */
Word16 freq_prev[MA_NP][M],
/* (i/o) : previous LPS for quantization        */
Word16 *seed, /* (i/o) : random generator seed                */
Word32 L_exc_err[4]) {

	Word16 i;

	Word16 curAcf[MP1];
	Word16 bid[M], zero[MP1];
	Word16 curCoeff[MP1];
	Word16 lsp_new[M];
	Word16 *lpcCoeff;
	Word16 cur_igain;
	Word16 energyq, temp;

	/* Update Ener and sh_ener */
	for (i = NB_GAIN - 1; i >= 1; i--) {
		dtxc->ener[i] = dtxc->ener[i - 1];
		dtxc->sh_ener[i] = dtxc->sh_ener[i - 1];
	}

	/* Compute current Acfs */
	Calc_sum_acf(dtxc->Acf, dtxc->sh_Acf, curAcf, &dtxc->sh_ener[0], NB_CURACF);

	/* Compute LPC coefficients and residual energy */
	if (curAcf[0] == 0) {
		dtxc->ener[0] = 0; /* should not happen */
	} else {
		Set_zero(zero, MP1);
		Levinson(dtxc, curAcf, zero, curCoeff, bid, &dtxc->ener[0]);
	}

	/* if first frame of silence => SID frame */
	if (pastVad != 0) {
		ana[0] = 2;
		dtxc->count_fr0 = 0;
		dtxc->nb_ener = 1;
		Qua_Sidgain(dtxc->ener, dtxc->sh_ener, dtxc->nb_ener, &energyq,
				&cur_igain);

	} else {
		dtxc->nb_ener = add(dtxc->nb_ener, 1);
		if (sub(dtxc->nb_ener, NB_GAIN) > 0)
			dtxc->nb_ener = NB_GAIN;
		Qua_Sidgain(dtxc->ener, dtxc->sh_ener, dtxc->nb_ener, &energyq,
				&cur_igain);

		/* Compute stationarity of current filter   */
		/* versus reference filter                  */
		if (Cmp_filt(dtxc->RCoeff, dtxc->sh_RCoeff, curAcf, dtxc->ener[0],
		FRAC_THRESH1) != 0) {
			dtxc->flag_chang = 1;
		}

		/* compare energy difference between current frame and last frame */
		temp = abs_s(sub(dtxc->prev_energy, energyq));
		temp = sub(temp, 2);
		if (temp > 0)
			dtxc->flag_chang = 1;

		dtxc->count_fr0 = add(dtxc->count_fr0, 1);
		if (sub(dtxc->count_fr0, FR_SID_MIN) < 0) {
			ana[0] = 0; /* no transmission */
		} else {
			if (dtxc->flag_chang != 0) {
				ana[0] = 2; /* transmit SID frame */
			} else {
				ana[0] = 0;
			}

			dtxc->count_fr0 = FR_SID_MIN; /* to avoid overflow */
		}
	}

	if (sub(ana[0], 2) == 0) {

		/* Reset frame count and change flag */
		dtxc->count_fr0 = 0;
		dtxc->flag_chang = 0;

		/* Compute past average filter */
		Calc_pastfilt(dtxc);
		Calc_RCoeff(dtxc->pastCoeff, dtxc->RCoeff, &dtxc->sh_RCoeff);

		/* Compute stationarity of current filter   */
		/* versus past average filter               */

		/* if stationary */
		/* transmit average filter => new ref. filter */
		if (Cmp_filt(dtxc->RCoeff, dtxc->sh_RCoeff, curAcf, dtxc->ener[0],
		FRAC_THRESH2) == 0) {
			lpcCoeff = dtxc->pastCoeff;
		}

		/* else */
		/* transmit current filter => new ref. filter */
		else {
			lpcCoeff = curCoeff;
			Calc_RCoeff(curCoeff, dtxc->RCoeff, &dtxc->sh_RCoeff);
		}

		/* Compute SID frame codes */

		Az_lsp(lpcCoeff, lsp_new, lsp_old_q); /* From A(z) to lsp */

		/* LSP quantization */
		lsfq_noise(lsp_new, dtxc->lspSid_q, freq_prev, &ana[1]);

		dtxc->prev_energy = energyq;
		ana[4] = cur_igain;
		dtxc->sid_gain = tab_Sidgain[cur_igain];

	} /* end of SID frame case */

	/* Compute new excitation */
	if (pastVad != 0) {
		dtxc->cur_gain = dtxc->sid_gain;
	} else {
		dtxc->cur_gain = mult_r(dtxc->cur_gain, A_GAIN0);
		dtxc->cur_gain = add(dtxc->cur_gain, mult_r(dtxc->sid_gain, A_GAIN1));
	}

	Calc_exc_rand(dtxc->cur_gain, exc, seed, L_exc_err);

	Int_qlpc(lsp_old_q, dtxc->lspSid_q, Aq);
	for (i = 0; i < M; i++) {
		lsp_old_q[i] = dtxc->lspSid_q[i];
	}

	/* Update sumAcf if fr_cur = 0 */
	if (dtxc->fr_cur == 0) {
		Update_sumAcf(dtxc);
	}

	return;
}

/*-----------------------------------------------------------*
 * procedure Update_cng:                                     *
 *           ~~~~~~~~~~                                      *
 *   Updates autocorrelation arrays                          *
 *   used for DTX/CNG                                        *
 *   If Vad=1 : updating of array sumAcf                     *
 *-----------------------------------------------------------*/
void Update_cng(struct dtx_ctx * restrict dtxc, Word16 * restrict r_h, /* (i) :   MSB of frame autocorrelation        */
Word16 exp_r, /* (i) :   scaling factor associated           */
Word16 Vad /* (i) :   current Vad decision                */
) {
	Word16 i;
	Word16 *ptr1, *ptr2;

	/* Update Acf and shAcf */
	ptr1 = dtxc->Acf + SIZ_ACF - 1;
	ptr2 = ptr1 - MP1;
	for (i = 0; i < (SIZ_ACF - MP1); i++) {
		*ptr1-- = *ptr2--;
	}
	for (i = NB_CURACF - 1; i >= 1; i--) {
		dtxc->sh_Acf[i] = dtxc->sh_Acf[i - 1];
	}

	/* Save current Acf */
	dtxc->sh_Acf[0] = negate(add(16, exp_r));
	for (i = 0; i < MP1; i++) {
		dtxc->Acf[i] = r_h[i];
	}

	dtxc->fr_cur = add(dtxc->fr_cur, 1);
	if (sub(dtxc->fr_cur, NB_CURACF) == 0) {
		dtxc->fr_cur = 0;
		if (Vad != 0) {
			Update_sumAcf(dtxc);
		}
	}

	return;
}

/*-----------------------------------------------------------*
 *         Local procedures                                  *
 *         ~~~~~~~~~~~~~~~~                                  *
 *-----------------------------------------------------------*/

/* Compute scaled autocorr of LPC coefficients used for Itakura distance */
/*************************************************************************/
static void Calc_RCoeff(Word16 *Coeff, Word16 *RCoeff, Word16 *sh_RCoeff) {
	Word16 i, j;
	Word16 sh1;
	Word32 L_acc;

	/* RCoeff[0] = SUM(j=0->M) Coeff[j] ** 2 */
	L_acc = 0L;
	for (j = 0; j <= M; j++) {
		L_acc = L_mac(L_acc, Coeff[j], Coeff[j]);
	}

	/* Compute exponent RCoeff */
	sh1 = norm_l(L_acc);
	L_acc = L_shl(L_acc, sh1);
	RCoeff[0] = g_round(L_acc);

	/* RCoeff[i] = SUM(j=0->M-i) Coeff[j] * Coeff[j+i] */
	for (i = 1; i <= M; i++) {
		L_acc = 0L;
		for (j = 0; j <= M - i; j++) {
			L_acc = L_mac(L_acc, Coeff[j], Coeff[j + i]);
		}
		L_acc = L_shl(L_acc, sh1);
		RCoeff[i] = g_round(L_acc);
	}
	*sh_RCoeff = sh1;
	return;
}

/* Compute Itakura distance and compare to threshold */
/*****************************************************/
static Word16 Cmp_filt(Word16 *RCoeff, Word16 sh_RCoeff, Word16 *acf,
		Word16 alpha, Word16 FracThresh) {
	Word32 L_temp0, L_temp1;
	Word16 temp1, temp2, sh[2], ind;
	Word16 i;
	Word16 diff, flag;

	sh[0] = 0;
	sh[1] = 0;
	ind = 1;
	flag = 0;
	do {
		clr_sat();
		temp1 = shr(RCoeff[0], sh[0]);
		temp2 = shr(acf[0], sh[1]);
		L_temp0 = L_shr(L_mult(temp1, temp2), 1);
		for (i = 1; i <= M; i++) {
			temp1 = shr(RCoeff[i], sh[0]);
			temp2 = shr(acf[i], sh[1]);
			L_temp0 = L_mac(L_temp0, temp1, temp2);
		}
		if (get_sat()) {
			sh[(int) ind] = add(sh[(int) ind], 1);
			ind = sub(1, ind);
		} else
			flag = 1;
	} while (flag == 0);

	temp1 = mult_r(alpha, FracThresh);
	L_temp1 = L_add(L_deposit_l(temp1), L_deposit_l(alpha));
	temp1 = add(sh_RCoeff, 9); /* 9 = Lpc_justif. * 2 - 16 + 1 */
	temp2 = add(sh[0], sh[1]);
	temp1 = sub(temp1, temp2);
	L_temp1 = L_shl(L_temp1, temp1);

	L_temp0 = L_sub(L_temp0, L_temp1);
	if (L_temp0 > 0L)
		diff = 1;
	else
		diff = 0;

	return (diff);
}

/* Compute past average filter */
/*******************************/
static void Calc_pastfilt(struct dtx_ctx *dtxc) {
	Word16 i;
	Word16 s_sumAcf[MP1];
	Word16 bid[M], zero[MP1];
	Word16 temp;

	Calc_sum_acf(dtxc->sumAcf, dtxc->sh_sumAcf, s_sumAcf, &temp, NB_SUMACF);

	if (s_sumAcf[0] == 0L) {
		dtxc->pastCoeff[0] = 4096;
		for (i = 1; i <= M; i++)
			dtxc->pastCoeff[i] = 0;
		return;
	}

	Set_zero(zero, MP1);
	Levinson(dtxc, s_sumAcf, zero, dtxc->pastCoeff, bid, &temp);
	return;
}

/* Update sumAcf */
/*****************/
static void Update_sumAcf(struct dtx_ctx *dtxc) {
	Word16 *ptr1, *ptr2;
	Word16 i;

	/*** Move sumAcf ***/
	ptr1 = dtxc->sumAcf + SIZ_SUMACF - 1;
	ptr2 = ptr1 - MP1;
	for (i = 0; i < (SIZ_SUMACF - MP1); i++) {
		*ptr1-- = *ptr2--;
	}
	for (i = NB_SUMACF - 1; i >= 1; i--) {
		dtxc->sh_sumAcf[i] = dtxc->sh_sumAcf[i - 1];
	}

	/* Compute new sumAcf */
	Calc_sum_acf(dtxc->Acf, dtxc->sh_Acf, dtxc->sumAcf, dtxc->sh_sumAcf,
	NB_CURACF);
	return;
}

/* Compute sum of acfs (curAcf, sumAcf or s_sumAcf) */
/****************************************************/
static void Calc_sum_acf(Word16 *acf, Word16 *sh_acf, Word16 *sum,
		Word16 *sh_sum, Word16 nb) {

	Word16 *ptr1;
	Word32 L_temp, L_tab[MP1];
	Word16 sh0, temp;
	Word16 i, j;

	/* Compute sum = sum of nb acfs */
	/* Find sh_acf minimum */
	sh0 = sh_acf[0];
	for (i = 1; i < nb; i++) {
		if (sub(sh_acf[i], sh0) < 0)
			sh0 = sh_acf[i];
	}
	sh0 = add(sh0, 14); /* 2 bits of margin */

	for (j = 0; j < MP1; j++) {
		L_tab[j] = 0L;
	}
	ptr1 = acf;
	for (i = 0; i < nb; i++) {
		temp = sub(sh0, sh_acf[i]);
		for (j = 0; j < MP1; j++) {
			L_temp = L_deposit_l(*ptr1++);
			L_temp = L_shl(L_temp, temp); /* shift right if temp<0 */
			L_tab[j] = L_add(L_tab[j], L_temp);
		}
	}
	temp = norm_l(L_tab[0]);
	for (i = 0; i <= M; i++) {
		sum[i] = extract_h(L_shl(L_tab[i], temp));
	}
	temp = sub(temp, 16);
	*sh_sum = add(sh0, temp);
	return;
}
