/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder with Annex B    ANSI-C Source Code
 Version 1.3    Last modified: August 1997

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies,
 Rockwell International
 All rights reserved.
 */

#include "typedef.h"
#include "ld8a.h"
#include "basic_op.h"
#include "tab_ld8a.h"

/* local function */
static inline void Lsp_iqua_cs(struct lspdec_ctx *lspdc, Word16 prm[], /* (i)     : indexes of the selected LSP */
                               Word16 lsp_q[], /* (o) Q13 : Quantized LSP parameters    */
                               Word16 erase /* (i)     : frame erase information     */
                               );

static inline void Mcopy(const Word16 *from, Word16 *to)
{
    int i;
    for (i = 0; i < M; i++)
        to[i] = from[i];
}

/*----------------------------------------------------------------------------
 * Lsp_decw_reset -   set the previous LSP vectors
 *----------------------------------------------------------------------------
 */
void Lsp_decw_reset(struct lspdec_ctx *lspdc)
{
    Word16 i;

    for (i = 0; i < MA_NP; i++)
        Mcopy(&freq_prev_reset[0], &lspdc->freq_prev[i][0]);

    lspdc->prev_ma = 0;

    Mcopy(freq_prev_reset, lspdc->prev_lsp);
}

/*----------------------------------------------------------------------------
 * Lsp_iqua_cs -  LSP main quantization routine
 *----------------------------------------------------------------------------
 */
static inline void Lsp_iqua_cs(struct lspdec_ctx *lspdc, Word16 prm[], /* (i)     : indexes of the selected LSP */
                               Word16 lsp_q[], /* (o) Q13 : Quantized LSP parameters    */
                               Word16 erase /* (i)     : frame erase information     */
                               )
{
    Word16 mode_index;
    Word16 code0;
    Word16 code1;
    Word16 code2;
    Word16 buf[M]; /* Q13 */

    if (erase == 0)
    { /* Not frame erasure */
        mode_index = shr(prm[0], NC0_B) & (Word16) 1;
        code0 = prm[0] & (Word16) (NC0 - 1);
        code1 = shr(prm[1], NC1_B) & (Word16) (NC1 - 1);
        code2 = prm[1] & (Word16) (NC1 - 1);

        /* compose quantized LSP (lsp_q) from indexes */

        Lsp_get_quant(lspcb1, lspcb2, code0, code1, code2, fg[mode_index],
                      lspdc->freq_prev, lsp_q, fg_sum[mode_index]);

        /* save parameters to use in case of the frame erased situation */

        Mcopy(lsp_q, lspdc->prev_lsp);
        lspdc->prev_ma = mode_index;
    }
    else
    { /* Frame erased */
        /* use revious LSP */

        Mcopy(lspdc->prev_lsp, lsp_q);

        /* update freq_prev */

        Lsp_prev_extract(lspdc->prev_lsp, buf, fg[lspdc->prev_ma],
                         (const Word16 (*)[M]) lspdc->freq_prev,
                         fg_sum_inv[lspdc->prev_ma]);
        Lsp_prev_update(buf, lspdc->freq_prev);
    }

    return;
}

/*-------------------------------------------------------------------*
 * Function  D_lsp:                                                  *
 *           ~~~~~~                                                  *
 *-------------------------------------------------------------------*/

void D_lsp(struct lspdec_ctx *lspdc, Word16 prm[], /* (i)     : indexes of the selected LSP */
           Word16 lsp_q[], /* (o) Q15 : Quantized LSP parameters    */
           Word16 erase /* (i)     : frame erase information     */
           )
{
    Word16 lsf_q[M]; /* domain 0.0<= lsf_q <PI in Q13 */

    Lsp_iqua_cs(lspdc, prm, lsf_q, erase);

    /* Convert LSFs to LSPs */

    Lsf_lsp2(lsf_q, lsp_q);

    return;
}

void Get_decfreq_prev(struct lspdec_ctx *lspdc, Word16 x[MA_NP][M])
{
    Word16 i;

    for (i = 0; i < MA_NP; i++)
        Mcopy(&lspdc->freq_prev[i][0], &x[i][0]);
}

void Update_decfreq_prev(struct lspdec_ctx *lspdc, Word16 x[MA_NP][M])
{
    Word16 i;

    for (i = 0; i < MA_NP; i++)
        Mcopy(&x[i][0], &lspdc->freq_prev[i][0]);
}

/*
 update previous LSP parameter
 */
void Lsp_prev_update(const Word16 lsp_ele[M], /* (i)   Q13 : LSP vectors           */
Word16 freq_prev[MA_NP][M] /* (i/o) Q13 : previous LSP vectors  */
) {
    Word16 k;

    for (k = MA_NP - 1; k > 0; k--)
        Mcopy(freq_prev[k - 1], freq_prev[k]);

    Mcopy(lsp_ele, freq_prev[0]);
    return;
}

