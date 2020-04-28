/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder    ANSI-C Source Code
 Version 1.1    Last modified: September 1996
 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies
 All rights reserved.
 */

/*-------------------------------------------------------------------*
 * Function  Convolve:                                               *
 *           ~~~~~~~~~                                               *
 *-------------------------------------------------------------------*
 * Perform the convolution between two vectors x[] and h[] and       *
 * write the result in the vector y[].                               *
 * All vectors are of length N.                                      *
 *-------------------------------------------------------------------*/

#include "typedef.h"
#include "basic_op.h"
#include "ld8a.h"

/*-----------------------------------------------------*
 * procedure Syn_filt:                                 *
 *           ~~~~~~~~                                  *
 * Do the synthesis filtering 1/A(z).                  *
 *-----------------------------------------------------*/
/* iterate L_SUBFR times */
void Syn_filt(Word16 a[], /* (i) Q12 : a[m+1] prediction coefficients   (m=10)  */
Word16 x[], /* (i)     : input signal                             */
Word16 y[restrict], /* (o)     : output signal                            */
Word16 mem[restrict], /* (i/o)   : memory associated with this filtering.   */
Word16 update /* (i)     : 0=no update, 1=update of memory.         */
) {
    Word16 i;
    Word32 s;
    Word16 tmp[100]; /* This is usually done by memory allocation (lg+M) */
    Word16 *yy;
    Word16 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    Word16 yy1, yy2, yy3, yy4, yy5, yy6, yy7, yy8, yy9, yy10;

    /* Copy mem[] to yy[] */

    yy = tmp;

    for (i = 0; i < M; i++) {
        *yy++ = mem[i];
    }

    /* Do the filtering. */
    a0 = a[0];
    a1 = a[1];
    a2 = a[2];
    a3 = a[3];
    a4 = a[4];
    a5 = a[5];
    a6 = a[6];
    a7 = a[7];
    a8 = a[8];
    a9 = a[9];
    a10 = a[10];
    yy1 = yy[-1];
    yy2 = yy[-2];
    yy3 = yy[-3];
    yy4 = yy[-4];
    yy5 = yy[-5];
    yy6 = yy[-6];
    yy7 = yy[-7];
    yy8 = yy[-8];
    yy9 = yy[-9];
    yy10 = yy[-10];

    for (i = 0; i < L_SUBFR; i++) {
        s = L_mult(x[i], a0);

        s = L_msu(s, a10, yy10);
        yy10 = yy9;
        s = L_msu(s, a9, yy9);
        yy9 = yy8;
        s = L_msu(s, a8, yy8);
        yy8 = yy7;
        s = L_msu(s, a7, yy7);
        yy7 = yy6;
        s = L_msu(s, a6, yy6);
        yy6 = yy5;
        s = L_msu(s, a5, yy5);
        yy5 = yy4;
        s = L_msu(s, a4, yy4);
        yy4 = yy3;
        s = L_msu(s, a3, yy3);
        yy3 = yy2;
        s = L_msu(s, a2, yy2);
        yy2 = yy1;
        s = L_msu(s, a1, yy1);

        s = L_shl(s, 3);
        yy1 = round(s);
        *yy++ = yy1;
    }

    for (i = 0; i < L_SUBFR; i++) {
        y[i] = tmp[i + M];
    }

    /* Update of memory if update==1 */

    if (update != 0)
        for (i = 0; i < M; i++) {
            mem[i] = y[L_SUBFR - M + i];
        }

    return;
}

/* iterate L_H times */
void Syn_filt_for_post(Word16 a[], /* (i) Q12 : a[m+1] prediction coefficients   (m=10)  */
Word16 x[], /* (i)     : input signal                             */
Word16 y[restrict], /* (o)     : output signal                            */
Word16 mem[restrict], /* (i/o)   : memory associated with this filtering.   */
Word16 update /* (i)     : 0=no update, 1=update of memory.         */
) {
    Word16 i;
    Word32 s;
    Word16 tmp[100]; /* This is usually done by memory allocation (lg+M) */
    Word16 *yy;
    Word16 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    Word16 yy1, yy2, yy3, yy4, yy5, yy6, yy7, yy8, yy9, yy10;

    /* Copy mem[] to yy[] */

    yy = tmp;

    for (i = 0; i < M; i++) {
        *yy++ = mem[i];
    }

    /* Do the filtering. */
    a0 = a[0];
    a1 = a[1];
    a2 = a[2];
    a3 = a[3];
    a4 = a[4];
    a5 = a[5];
    a6 = a[6];
    a7 = a[7];
    a8 = a[8];
    a9 = a[9];
    a10 = a[10];
    yy1 = yy[-1];
    yy2 = yy[-2];
    yy3 = yy[-3];
    yy4 = yy[-4];
    yy5 = yy[-5];
    yy6 = yy[-6];
    yy7 = yy[-7];
    yy8 = yy[-8];
    yy9 = yy[-9];
    yy10 = yy[-10];

    for (i = 0; i < L_H; i++) {
        s = L_mult(x[i], a0);

        s = L_msu(s, a10, yy10);
        yy10 = yy9;
        s = L_msu(s, a9, yy9);
        yy9 = yy8;
        s = L_msu(s, a8, yy8);
        yy8 = yy7;
        s = L_msu(s, a7, yy7);
        yy7 = yy6;
        s = L_msu(s, a6, yy6);
        yy6 = yy5;
        s = L_msu(s, a5, yy5);
        yy5 = yy4;
        s = L_msu(s, a4, yy4);
        yy4 = yy3;
        s = L_msu(s, a3, yy3);
        yy3 = yy2;
        s = L_msu(s, a2, yy2);
        yy2 = yy1;
        s = L_msu(s, a1, yy1);

        s = L_shl(s, 3);
        yy1 = round(s);
        *yy++ = yy1;
    }

    for (i = 0; i < L_H; i++) {
        y[i] = tmp[i + M];
    }

    /* Update of memory if update==1 */

    if (update != 0)
        for (i = 0; i < M; i++) {
            mem[i] = y[L_H - M + i];
        }

    return;
}

/*-----------------------------------------------------------------------*
 * procedure Residu:                                                     *
 *           ~~~~~~                                                      *
 * Compute the LPC residual  by filtering the input speech through A(z)  *
 *-----------------------------------------------------------------------*/

void Residu(Word16 a[restrict], /* (i) Q12 : prediction coefficients                     */
Word16 x[restrict], /* (i)     : speech (values x[-m..-1] are needed         */
Word16 y[restrict] /* (o)     : residual signal                             */
) {
	Word16 i, j;
	Word32 s;

	for (i = 0; i < L_SUBFR; i++) {
		s = L_mult(x[i], a[0]);
		for (j = 1; j <= M; j++)
			s = L_mac(s, a[j], x[i - j]);

		s = L_shl(s, 3);
		y[i] = g_round(s);
	}
	return;
}
