/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder    ANSI-C Source Code
 Version 1.1    Last modified: September 1996

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies
 All rights reserved.
 */

/*------------------------------------------------------------------------*
 * Function Pre_Process()                                                 *
 *                                                                        *
 * Preprocessing of input speech.                                         *
 *   - 2nd order high pass filter with cut off frequency at 140 Hz.       *
 *   - Divide input by two.                                               *
 *-----------------------------------------------------------------------*/

#include "typedef.h"
#include "basic_op.h"
#include "oper_32b.h"

#include "ld8a.h"
#include "tab_ld8a.h"

/*------------------------------------------------------------------------*
 * 2nd order high pass filter with cut off frequency at 140 Hz.           *
 * Designed with SPPACK efi command -40 dB att, 0.25 ri.                  *
 *                                                                        *
 * Algorithm:                                                             *
 *                                                                        *
 *  y[i] = b[0]*x[i]/2 + b[1]*x[i-1]/2 + b[2]*x[i-2]/2                    *
 *                     + a[1]*y[i-1]   + a[2]*y[i-2];                     *
 *                                                                        *
 *     b[3] = {0.92727435E+00, -0.18544941E+01, 0.92727435E+00};          *
 *     a[3] = {0.10000000E+01, 0.19059465E+01, -0.91140240E+00};          *
 *                                                                        *
 *  Input are divided by two in the filtering process.                    *
 *-----------------------------------------------------------------------*/

void Init_Pre_Process(struct prepost_process_ctx *preprc) {
	preprc->y2_hi = 0;
	preprc->y2_lo = 0;
	preprc->y1_hi = 0;
	preprc->y1_lo = 0;
	preprc->x0 = 0;
	preprc->x1 = 0;
}

void Pre_Process(struct prepost_process_ctx *preprc, Word16 signal[])
{
    Word16 i, x2;
    Word32 L_tmp, L_tmp1;
    Word16 r_x0, r_x1;
    Word16 r_y1_hi, r_y1_lo, r_y2_hi, r_y2_lo;
    Word16 a140_1, a140_2;
    Word16 b140_0, b140_1, b140_2;

    r_x1 = preprc->x1;
    r_x0 = preprc->x0;
    r_y1_hi = preprc->y1_hi;
    r_y1_lo = preprc->y1_lo;
    r_y2_hi = preprc->y2_hi;
    r_y2_lo = preprc->y2_lo;
    a140_1 = a140[1];
    a140_2 = a140[2];
    b140_0 = b140[0];
    b140_1 = b140[1];
    b140_2 = b140[2];

    for (i = 0; i < L_FRAME; i++)
    {
        x2 = r_x1;
        r_x1 = r_x0;
        r_x0 = signal[i];

        /*  y[i] = b[0]*x[i]/2 + b[1]*x[i-1]/2 + b140[2]*x[i-2]/2  */
        /*                     + a[1]*y[i-1] + a[2] * y[i-2];      */

        L_tmp = L_mult(r_y1_hi, a140_1);
        L_tmp = L_mac(L_tmp, mult(r_y1_lo, a140_1), 1);
        L_tmp1 = L_mult(r_y2_hi, a140_2);
        L_tmp1 = L_mac(L_tmp1, mult(r_y2_lo, a140_2), 1);
        L_tmp = L_add(L_tmp, L_tmp1);

        L_tmp = L_mac(L_tmp, r_x0, b140_0);
        L_tmp = L_mac(L_tmp, r_x1, b140_1);
        L_tmp = L_mac(L_tmp, x2, b140_2);
        L_tmp = L_shl(L_tmp, 3); /* Q28 --> Q31 (Q12 --> Q15) */
        signal[i] = g_round(L_tmp);

        r_y2_hi = r_y1_hi;
        r_y2_lo = r_y1_lo;
        L_Extract(L_tmp, &r_y1_hi, &r_y1_lo);
    }

    preprc->y1_hi = r_y1_hi;
    preprc->y1_lo = r_y1_lo;
    preprc->y2_hi = r_y2_hi;
    preprc->y2_lo = r_y2_lo;
    preprc->x1 = r_x1;
    preprc->x0 = r_x0;

    return;
}
