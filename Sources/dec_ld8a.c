/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder with Annex B    ANSI-C Source Code
 Version 1.5    Last modified: October 2006

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies,
 Rockwell International
 All rights reserved.
 */

/*-----------------------------------------------------------------*
 *   Functions Init_Decod_ld8a  and Decod_ld8a                     *
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "typedef.h"
#include "basic_op.h"
#include "ld8a.h"

#include "sid.h"

/*---------------------------------------------------------------*
 *   Decoder constant parameters (defined in "ld8a.h")           *
 *---------------------------------------------------------------*
 *   L_FRAME     : Frame size.                                   *
 *   L_SUBFR     : Sub-frame size.                               *
 *   M           : LPC order.                                    *
 *   MP1         : LPC order+1                                   *
 *   PIT_MIN     : Minimum pitch lag.                            *
 *   PIT_MAX     : Maximum pitch lag.                            *
 *   L_INTERPOL  : Length of filter for interpolation            *
 *   PRM_SIZE    : Size of vector containing analysis parameters *
 *---------------------------------------------------------------*/

/* Lsp (Line spectral pairs) */

static const Word16 lsp_old_c[M] = { 30000, 26000, 21000, 15000, 8000, 0, -8000,
		-15000, -21000, -26000 };

/* Gain predictor, Past quantized energies = -14.0 in Q10 */
static Word16 past_qua_en_c[4] = { -14336, -14336, -14336, -14336 };

/*-----------------------------------------------------------------*
 *   Function Init_Decod_ld8a                                      *
 *            ~~~~~~~~~~~~~~~                                      *
 *                                                                 *
 *   ->Initialization of variables for the decoder section.        *
 *                                                                 *
 *-----------------------------------------------------------------*/

void Init_Decod_ld8a(struct dec_ld8a_ctx *decc) {

	/* Initialize static pointer */

	decc->exc = decc->old_exc + PIT_MAX + L_INTERPOL;

	/* Static vectors to zero */

	Set_zero(decc->old_exc, PIT_MAX + L_INTERPOL);
	Set_zero(decc->mem_syn, M);

	decc->sharp = SHARPMIN;
	decc->old_T0 = 60;
	decc->gain_code = 0;
	decc->gain_pitch = 0;

	Copy(lsp_old_c, decc->lsp_old, M);
	Copy(past_qua_en_c, decc->past_qua_en, 4);

	Lsp_decw_reset(&decc->lspdc);

	/* for G.729B */
	decc->seed_fer = 21845;
	decc->past_ftyp = 1;
	decc->seed = INIT_SEED;
	decc->sid_sav = 0;
	decc->sh_sid_sav = 1;

	return;
}

/*-----------------------------------------------------------------*
 *   Function Decod_ld8a                                           *
 *           ~~~~~~~~~~                                            *
 *   ->Main decoder routine.                                       *
 *                                                                 *
 *-----------------------------------------------------------------*/

void Decod_ld8a(struct dec_ld8a_ctx *decc, Word16 parm[], /* (i)   : vector of synthesis parameters
 parm[0] = bad frame indicator (bfi)  */
Word16 synth[], /* (o)   : synthesis speech                     */
Word16 A_t[], /* (o)   : decoded LP filter in 2 subframes     */
Word16 *T2, /* (o)   : decoded pitch lag in 2 subframes     */
Word16 *Vad /* (o)   : frame type                           */
) {
	Word16 *Az; /* Pointer on A_t   */
	Word16 lsp_new[M]; /* LSPs             */
	Word16 code[L_SUBFR]; /* ACELP codevector */

	/* Scalars */

	Word16 i, j, i_subfr;
	Word16 T0, T0_frac, index;
	Word16 bfi;
	Word32 L_temp;

	Word16 bad_pitch; /* bad pitch indicator */

	/* for G.729B */
	Word16 ftyp;
	Word16 lsfq_mem[MA_NP][M];

	/* Test bad frame indicator (bfi) */

	bfi = *parm++;
	/* for G.729B */
	ftyp = *parm;

	if (bfi == 1) {
		if (decc->past_ftyp == 1) {
			ftyp = 1;
			parm[4] = 1; /* G.729 maintenance */
		} else
			ftyp = 0;
		*parm = ftyp; /* modification introduced in version V1.3 */
	}

	*Vad = ftyp;

	/* Processing non active frames (SID & not transmitted) */
	if (ftyp != 1) {

		Get_decfreq_prev(&decc->lspdc, lsfq_mem);
		Dec_cng(&decc->decsidc, decc->past_ftyp, decc->sid_sav,
				decc->sh_sid_sav, parm, decc->exc, decc->lsp_old, A_t,
				&decc->seed, lsfq_mem);
		Update_decfreq_prev(&decc->lspdc, lsfq_mem);

		Az = A_t;
		for (i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR) {
			clr_sat();
			Syn_filt(Az, &decc->exc[i_subfr], &synth[i_subfr],
					decc->mem_syn, 0);
			if (get_sat()) {
				/* In case of overflow in the synthesis          */
				/* -> Scale down vector exc[] and redo synthesis */

				for (i = 0; i < PIT_MAX + L_INTERPOL + L_FRAME; i++)
					decc->old_exc[i] = shr(decc->old_exc[i], 2);

				Syn_filt(Az, &decc->exc[i_subfr], &synth[i_subfr],
						decc->mem_syn, 1);
			} else
				Copy(&synth[i_subfr + L_SUBFR - M], decc->mem_syn, M);

			Az += MP1;

			*T2++ = decc->old_T0;
		}
		decc->sharp = SHARPMIN;

	}
	/* Processing active frame */
	else {

		decc->seed = INIT_SEED;
		parm++;

		/* Decode the LSPs */

		D_lsp(&decc->lspdc, parm, lsp_new, bfi);
		parm += 2;

		/* Interpolation of LPC for the 2 subframes */

		Int_qlpc(decc->lsp_old, lsp_new, A_t);

		/* update the LSFs for the next frame */

		Copy(lsp_new, decc->lsp_old, M);

		/*------------------------------------------------------------------------*
		 *          Loop for every subframe in the analysis frame                 *
		 *------------------------------------------------------------------------*
		 * The subframe size is L_SUBFR and the loop is repeated L_FRAME/L_SUBFR  *
		 *  times                                                                 *
		 *     - decode the pitch delay                                           *
		 *     - decode algebraic code                                            *
		 *     - decode pitch and codebook gains                                  *
		 *     - find the excitation and compute synthesis speech                 *
		 *------------------------------------------------------------------------*/

		Az = A_t; /* pointer to interpolated LPC parameters */

		for (i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR)
		{

			index = *parm++; /* pitch index */

			if (i_subfr == 0) {
				i = *parm++; /* get parity check result */
				bad_pitch = add(bfi, i);
				if (bad_pitch == 0) {
					Dec_lag3(index, PIT_MIN, PIT_MAX, i_subfr, &T0, &T0_frac);
					decc->old_T0 = T0;
				} else /* Bad frame, or parity error */
				{
					T0 = decc->old_T0;
					T0_frac = 0;
					decc->old_T0 = add(decc->old_T0, 1);
					if (sub(decc->old_T0, PIT_MAX) > 0) {
						decc->old_T0 = PIT_MAX;
					}
				}
			} else /* second subframe */
			{
				if (bfi == 0) {
					Dec_lag3(index, PIT_MIN, PIT_MAX, i_subfr, &T0, &T0_frac);
					decc->old_T0 = T0;
				} else {
					T0 = decc->old_T0;
					T0_frac = 0;
					decc->old_T0 = add(decc->old_T0, 1);
					if (sub(decc->old_T0, PIT_MAX) > 0) {
						decc->old_T0 = PIT_MAX;
					}
				}
			}
			*T2++ = T0;

			/*-------------------------------------------------*
			 * - Find the adaptive codebook vector.            *
			 *-------------------------------------------------*/

			Pred_lt_3(&decc->exc[i_subfr], T0, T0_frac);

			/*-------------------------------------------------------*
			 * - Decode innovative codebook.                         *
			 * - Add the fixed-gain pitch contribution to code[].    *
			 *-------------------------------------------------------*/

			if (bfi != 0) /* Bad frame */
			{

				parm[0] = Random(&decc->seed_fer) & (Word16) 0x1fff; /* 13 bits random */
				parm[1] = Random(&decc->seed_fer) & (Word16) 0x000f; /*  4 bits random */
			}

			Decod_ACELP(parm[1], parm[0], code);
			parm += 2;

			j = shl(decc->sharp, 1); /* From Q14 to Q15 */
			if (sub(T0, L_SUBFR) < 0) {
				for (i = T0; i < L_SUBFR; i++) {
					code[i] = add(code[i], mult(code[i - T0], j));
				}
			}

			/*-------------------------------------------------*
			 * - Decode pitch and codebook gains.              *
			 *-------------------------------------------------*/

			index = *parm++; /* index of energy VQ */

			Dec_gain(index, code, bfi, &decc->gain_pitch,
					&decc->gain_code, decc->past_qua_en);

			/*-------------------------------------------------------------*
			 * - Update pitch sharpening "sharp" with quantized gain_pitch *
			 *-------------------------------------------------------------*/

			decc->sharp = decc->gain_pitch;
			if (sub(decc->sharp, SHARPMAX) > 0) {
				decc->sharp = SHARPMAX;
			}
			if (sub(decc->sharp, SHARPMIN) < 0) {
				decc->sharp = SHARPMIN;
			}

			/*-------------------------------------------------------*
			 * - Find the total excitation.                          *
			 * - Find synthesis speech corresponding to exc[].       *
			 *-------------------------------------------------------*/

			for (i = 0; i < L_SUBFR; i++) {
				/* exc[i] = gain_pitch*exc[i] + gain_code*code[i]; */
				/* exc[i]  in Q0   gain_pitch in Q14               */
				/* code[i] in Q13  gain_codeode in Q1              */

				L_temp = L_mult(decc->exc[i + i_subfr], decc->gain_pitch);
				L_temp = L_mac(L_temp, code[i], decc->gain_code);
				L_temp = L_shl(L_temp, 1);
				decc->exc[i + i_subfr] = g_round(L_temp);
			}

			clr_sat();
			Syn_filt(Az, &decc->exc[i_subfr], &synth[i_subfr],
					decc->mem_syn, 0);
			if (get_sat()) {
				/* In case of overflow in the synthesis          */
				/* -> Scale down vector exc[] and redo synthesis */

				for (i = 0; i < PIT_MAX + L_INTERPOL + L_FRAME; i++)
					decc->old_exc[i] = shr(decc->old_exc[i], 2);

				Syn_filt(Az, &decc->exc[i_subfr], &synth[i_subfr],
						decc->mem_syn, 1);
			} else
				Copy(&synth[i_subfr + L_SUBFR - M], decc->mem_syn, M);

			Az += MP1; /* interpolated LPC parameters for next subframe */
		}
	}

	/*------------*
	 *  For G729b
	 *-----------*/
	if (bfi == 0) {
		L_temp = 0L;
		for (i = 0; i < L_FRAME; i++) {
			L_temp = L_mac(L_temp, decc->exc[i], decc->exc[i]);
		} /* may overflow => last level of SID quantizer */
		decc->sh_sid_sav = norm_l(L_temp);
		decc->sid_sav = g_round(L_shl(L_temp, decc->sh_sid_sav));
		decc->sh_sid_sav = sub(16, decc->sh_sid_sav);
	}

	/*--------------------------------------------------*
	 * Update signal for next frame.                    *
	 * -> shift to the left by L_FRAME  exc[]           *
	 *--------------------------------------------------*/

	Copy(&decc->old_exc[L_FRAME], &decc->old_exc[0], PIT_MAX + L_INTERPOL);

	/* for G729b */
	decc->past_ftyp = ftyp;

	return;
}
