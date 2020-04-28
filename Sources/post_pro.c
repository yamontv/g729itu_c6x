/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder    ANSI-C Source Code
 Version 1.1    Last modified: September 1996

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies
 All rights reserved.
 */

/*------------------------------------------------------------------------*
 * Function Post_Process()                                                *
 *                                                                        *
 * Post-processing of output speech.                                      *
 *   - 2nd order high pass filter with cut off frequency at 100 Hz.       *
 *   - Multiplication by two of output speech with saturation.            *
 *-----------------------------------------------------------------------*/

#include "typedef.h"
#include "basic_op.h"
#include "oper_32b.h"

#include "ld8a.h"
#include "tab_ld8a.h"

/*------------------------------------------------------------------------*
 * 2nd order high pass filter with cut off frequency at 100 Hz.           *
 * Designed with SPPACK efi command -40 dB att, 0.25 ri.                  *
 *                                                                        *
 * Algorithm:                                                             *
 *                                                                        *
 *  y[i] = b[0]*x[i]   + b[1]*x[i-1]   + b[2]*x[i-2]                      *
 *                     + a[1]*y[i-1]   + a[2]*y[i-2];                     *
 *                                                                        *
 *     b[3] = {0.93980581E+00, -0.18795834E+01, 0.93980581E+00};          *
 *     a[3] = {0.10000000E+01, 0.19330735E+01, -0.93589199E+00};          *
 *-----------------------------------------------------------------------*/

/* Initialization of static values */

void Init_Post_Process(struct prepost_process_ctx *postprc) {
	postprc->y2_hi = 0;
	postprc->y2_lo = 0;
	postprc->y1_hi = 0;
	postprc->y1_lo = 0;
	postprc->x0 = 0;
	postprc->x1 = 0;
}

void Post_Process(struct prepost_process_ctx *postprc, Word16 signal[])
{
    Word16 i, x2;
    Word32 L_tmp, L_tmp1;
    Word16 r_x0, r_x1;
    Word16 r_y1_hi, r_y1_lo, r_y2_hi, r_y2_lo;
    Word16 a100_1, a100_2;
    Word16 b100_0, b100_1, b100_2;

    r_x1 = postprc->x1;
    r_x0 = postprc->x0;
    r_y1_hi = postprc->y1_hi;
    r_y1_lo = postprc->y1_lo;
    r_y2_hi = postprc->y2_hi;
    r_y2_lo = postprc->y2_lo;
    a100_1 = a100[1];
    a100_2 = a100[2];
    b100_0 = b100[0];
    b100_1 = b100[1];
    b100_2 = b100[2];

    for (i = 0; i < L_FRAME; i++)
    {
        x2 = r_x1;
        r_x1 = r_x0;
        r_x0 = signal[i];

        /*  y[i] = b[0]*x[i]   + b[1]*x[i-1]   + b[2]*x[i-2]    */
        /*                     + a[1]*y[i-1] + a[2] * y[i-2];      */

        L_tmp = L_mult(r_y1_hi, a100_1);
        L_tmp = L_mac(L_tmp, mult(r_y1_lo, a100_1), 1);
        L_tmp1 = L_mult(r_y2_hi, a100_2);
        L_tmp1 = L_mac(L_tmp1, mult(r_y2_lo, a100_2), 1);
        L_tmp = L_add(L_tmp, L_tmp1);

        L_tmp = L_mac(L_tmp, r_x0, b100_0);
        L_tmp = L_mac(L_tmp, r_x1, b100_1);
        L_tmp = L_mac(L_tmp, x2, b100_2);
        L_tmp = L_shl(L_tmp, 2); /* Q29 --> Q31 (Q13 --> Q15) */

        /* Multiplication by two of output speech with saturation. */
        signal[i] = g_round(L_shl(L_tmp, 1));

        r_y2_hi = r_y1_hi;
        r_y2_lo = r_y1_lo;
        L_Extract(L_tmp, &r_y1_hi, &r_y1_lo);
    }
    postprc->y1_hi = r_y1_hi;
    postprc->y1_lo = r_y1_lo;
    postprc->y2_hi = r_y2_hi;
    postprc->y2_lo = r_y2_lo;
    postprc->x1 = r_x1;
    postprc->x0 = r_x0;

    return;
}
