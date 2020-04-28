/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Speech Coder with Annex B    ANSI-C Source Code
 Version 1.3    Last modified: August 1997

 Copyright (c) 1996,
 AT&T, France Telecom, NTT, Universite de Sherbrooke, Lucent Technologies,
 Rockwell International
 All rights reserved.
 */

/*---------------------------------------------------------------*
 * LD8A.H                                                        *
 * ~~~~~~                                                        *
 * Function prototypes and constants use for G.729A 8kb/s coder. *
 *                                                               *
 *---------------------------------------------------------------*/

/*--------------------------------------------------------------------------*
 *       Codec constant parameters (coder, decoder, and postfilter)         *
 *--------------------------------------------------------------------------*/

#define  L_TOTAL      240     /* Total size of speech buffer.               */
#define  L_WINDOW     240     /* Window size in LP analysis.                */
#define  L_NEXT       40      /* Lookahead in LP analysis.                  */
#define  L_FRAME      80      /* Frame size.                                */
#define  L_SUBFR      40      /* Subframe size.                             */
#define  M            10      /* Order of LP filter.                        */
#define  MP1          (M+1)   /* Order of LP filter + 1                     */
#define  PIT_MIN      20      /* Minimum pitch lag.                         */
#define  PIT_MAX      143     /* Maximum pitch lag.                         */
#define  L_INTERPOL   (10+1)  /* Length of filter for interpolation.        */
#define  GAMMA1       24576   /* Bandwitdh factor = 0.75   in Q15           */

#define  PRM_SIZE     11      /* Size of vector of analysis parameters.     */
#define  SERIAL_SIZE  (80+2)  /* bfi+ number of speech bits                 */

#define SHARPMAX  13017   /* Maximum value of pitch sharpening     0.8  Q14 */
#define SHARPMIN  3277    /* Minimum value of pitch sharpening     0.2  Q14 */

/*--------------------------------------------------------------------------*
 *       LSP constant parameters                                            *
 *--------------------------------------------------------------------------*/

#define   NC            5      /*  NC = M/2 */
#define   MA_NP         4      /* MA prediction order for LSP */
#define   MODE          2      /* number of modes for MA prediction */
#define   NC0_B         7      /* number of first stage bits */
#define   NC1_B         5      /* number of second stage bits */
#define   NC0           (1<<NC0_B)
#define   NC1           (1<<NC1_B)

#define   L_LIMIT          40   /* Q13:0.005 */
#define   M_LIMIT       25681   /* Q13:3.135 */

#define   GAP1          10     /* Q13 */
#define   GAP2          5      /* Q13 */
#define   GAP3          321    /* Q13 */
#define GRID_POINTS     50

#define PI04      ((Word16)1029)        /* Q13  pi*0.04 */
#define PI92      ((Word16)23677)       /* Q13  pi*0.92 */
#define CONST10   ((Word16)10*(1<<11))  /* Q11  10.0 */
#define CONST12   ((Word16)19661)       /* Q14  1.2 */

/*-------------------------------*
 * Mathematic functions.         *
 *-------------------------------*/

Word32 Inv_sqrt( /* (o) Q30 : output value   (range: 0<=val<1)           */
Word32 L_x /* (i) Q0  : input value    (range: 0<=val<=7fffffff)   */
);

void Log2(Word32 L_x, /* (i) Q0 : input value                                 */
Word16 *exponent, /* (o) Q0 : Integer part of Log2.   (range: 0<=val<=30) */
Word16 *fraction /* (o) Q15: Fractionnal part of Log2. (range: 0<=val<1) */
);

Word32 Pow2( /* (o) Q0  : result       (range: 0<=val<=0x7fffffff) */
Word16 exponent, /* (i) Q0  : Integer part.      (range: 0<=val<=30)   */
Word16 fraction /* (i) Q15 : Fractionnal part.  (range: 0.0<=val<1.0) */
);

/*-------------------------------*
 * Pre and post-process.         *
 *-------------------------------*/
struct prepost_process_ctx {
	/* y[] values is keep in double precision      */
	Word16 y2_hi, y2_lo, y1_hi, y1_lo, x0, x1;
};

void Init_Pre_Process(struct prepost_process_ctx *preprc);

void Init_Post_Process(struct prepost_process_ctx *postprc);

void Pre_Process(struct prepost_process_ctx *preprc, Word16 signal[]);

void Post_Process(struct prepost_process_ctx *postprc, Word16 signal[]);

/*-------------------------------*
 * LPC analysis and filtering.   *
 *-------------------------------*/

void Autocorr(Word16 x[], /* (i)    : Input signal                      */
Word16 r_h[], /* (o)    : Autocorrelations  (msb)           */
Word16 r_l[], /* (o)    : Autocorrelations  (lsb)           */
Word16 *exp_R0);

void Lag_window(Word16 r_h[], /* (i/o)   : Autocorrelations  (msb)          */
Word16 r_l[] /* (i/o)   : Autocorrelations  (lsb)          */
);

void Az_lsp(Word16 a[], /* (i) Q12 : predictor coefficients              */
Word16 lsp[], /* (o) Q15 : line spectral pairs                 */
Word16 old_lsp[] /* (i)     : old lsp[] (in case not found 10 roots) */
);

void Lsp_Az(Word16 lsp[], /* (i) Q15 : line spectral frequencies            */
Word16 a[] /* (o) Q12 : predictor coefficients (order = 10)  */
);

void Lsp_lsf(Word16 lsp[], /* (i) Q15 : lsp[m] (range: -1<=val<1)                */
Word16 lsf[] /* (o) Q15 : lsf[m] normalized (range: 0.0<=val<=0.5) */
);

void Int_qlpc(Word16 lsp_old[], /* input : LSP vector of past frame              */
Word16 lsp_new[], /* input : LSP vector of present frame           */
Word16 Az[] /* output: interpolated Az() for the 2 subframes */
);

void Weight_Az(Word16 a[], /* (i) Q12 : a[m+1]  LPC coefficients             */
Word16 gamma, /* (i) Q15 : Spectral expansion factor.           */
Word16 ap[] /* (o) Q12 : Spectral expanded LPC coefficients   */
);

void Residu(Word16 a[], /* (i) Q12 : prediction coefficients                     */
Word16 x[], /* (i)     : speech (values x[-m..-1] are needed         */
Word16 y[] /* (o)     : residual signal                             */
);

void Syn_filt(Word16 a[], /* (i) Q12 : a[m+1] prediction coefficients   (m=10)  */
Word16 x[], /* (i)     : input signal                             */
Word16 y[], /* (o)     : output signal                            */
Word16 mem[], /* (i/o)   : memory associated with this filtering.   */
Word16 update /* (i)     : 0=no update, 1=update of memory.         */
);

void Syn_filt_for_post(Word16 a[], /* (i) Q12 : a[m+1] prediction coefficients   (m=10)  */
Word16 x[], /* (i)     : input signal                             */
Word16 y[], /* (o)     : output signal                            */
Word16 mem[], /* (i/o)   : memory associated with this filtering.   */
Word16 update /* (i)     : 0=no update, 1=update of memory.         */
);

/*--------------------------------------------------------------------------*
 *       LTP constant parameters                                            *
 *--------------------------------------------------------------------------*/

#define UP_SAMP         3
#define L_INTER10       10
#define FIR_SIZE_SYN    (UP_SAMP*L_INTER10+1)

/*-----------------------*
 * Pitch functions.      *
 *-----------------------*/

Word16 Pitch_ol_fast( /* output: open loop pitch lag                        */
Word16 signal[] /* input : signal used to compute the open loop pitch */
/*     signal[-pit_max] to signal[-1] should be known */
);

Word16 Pitch_fr3_fast(/* (o)     : pitch period.                          */
Word16 exc[], /* (i)     : excitation buffer                      */
Word16 xn[], /* (i)     : target vector                          */
Word16 h[], /* (i) Q12 : impulse response of filters.           */
Word16 t0_min, /* (i)     : minimum value in the searched range.   */
Word16 t0_max, /* (i)     : maximum value in the searched range.   */
Word16 i_subfr, /* (i)     : indicator for first subframe.          */
Word16 *pit_frac /* (o)     : chosen fraction.                       */
);

Word16 G_pitch( /* (o) Q14 : Gain of pitch lag saturated to 1.2       */
Word16 xn[], /* (i)     : Pitch target.                            */
Word16 y1[], /* (i)     : Filtered adaptive codebook.              */
Word16 g_coeff[] /* (i)     : Correlations need for gain quantization. */
);

Word16 Enc_lag3( /* output: Return index of encoding */
Word16 T0, /* input : Pitch delay              */
Word16 T0_frac, /* input : Fractional pitch delay   */
Word16 *T0_min, /* in/out: Minimum search delay     */
Word16 *T0_max, /* in/out: Maximum search delay     */
Word16 pit_min, /* input : Minimum pitch delay      */
Word16 pit_max, /* input : Maximum pitch delay      */
Word16 pit_flag /* input : Flag for 1st subframe    */
);

void Dec_lag3( /* output: return integer pitch lag       */
Word16 index, /* input : received pitch index           */
Word16 pit_min, /* input : minimum pitch lag              */
Word16 pit_max, /* input : maximum pitch lag              */
Word16 i_subfr, /* input : subframe flag                  */
Word16 *T0, /* output: integer part of pitch lag      */
Word16 *T0_frac /* output: fractional part of pitch lag   */
);

Word16 Interpol_3( /* (o)  : interpolated value  */
Word16 *x, /* (i)  : input vector        */
Word16 frac /* (i)  : fraction            */
);

void Pred_lt_3(Word16 exc[], /* in/out: excitation buffer */
Word16 T0, /* input : integer pitch lag */
Word16 frac /* input : fraction of lag   */
);

Word16 Parity_Pitch( /* output: parity bit (XOR of 6 MSB bits)    */
Word16 pitch_index /* input : index for which parity to compute */
);

Word16 Check_Parity_Pitch( /* output: 0 = no error, 1= error */
Word16 pitch_index, /* input : index of parameter     */
Word16 parity /* input : parity bit             */
);

void Cor_h_X(Word16 h[], /* (i) Q12 :Impulse response of filters      */
Word16 X[], /* (i)     :Target vector                    */
Word16 D[] /* (o)     :Correlations between h[] and D[] */
/*          Normalized to 13 bits            */
);

/*-----------------------*
 * Innovative codebook.  *
 *-----------------------*/

#define DIM_RR  616 /* size of correlation matrix                            */
#define NB_POS  8   /* Number of positions for each pulse                    */
#define STEP    5   /* Step betweem position of the same pulse.              */
#define MSIZE   64  /* Size of vectors for cross-correlation between 2 pulses*/

/* The following constants are Q15 fractions.
 These fractions is used to keep maximum precision on "alp" sum */

#define _1_2    (Word16)(16384)
#define _1_4    (Word16)( 8192)
#define _1_8    (Word16)( 4096)
#define _1_16   (Word16)( 2048)

Word16 ACELP_Code_A( /* (o)     :index of pulses positions    */
Word16 x[], /* (i)     :Target vector                */
Word16 h[], /* (i) Q12 :Inpulse response of filters  */
Word16 T0, /* (i)     :Pitch lag                    */
Word16 pitch_sharp, /* (i) Q14 :Last quantized pitch gain    */
Word16 code[], /* (o) Q13 :Innovative codebook          */
Word16 y[], /* (o) Q12 :Filtered innovative codebook */
Word16 *sign /* (o)     :Signs of 4 pulses            */
);

void Decod_ACELP(Word16 sign, /* (i)     : signs of 4 pulses.                       */
Word16 index, /* (i)     : Positions of the 4 pulses.               */
Word16 cod[] /* (o) Q13 : algebraic (fixed) codebook excitation    */
);

/*-------------------------------*
 * LSP VQ functions.             *
 *-------------------------------*/

void Lsf_lsp2(Word16 lsf[], /* (i) Q13 : lsf[m] (range: 0.0<=val<PI) */
Word16 lsp[] /* (o) Q15 : lsp[m] (range: -1<=val<1)   */
);

void Lsp_lsf2(Word16 lsp[], /* (i) Q15 : lsp[m] (range: -1<=val<1)   */
Word16 lsf[] /* (o) Q13 : lsf[m] (range: 0.0<=val<PI) */
);

void Qua_lsp(Word16 qua_freq_prev[MA_NP][M], Word16 lsp[], /* (i) Q15 : Unquantized LSP            */
Word16 lsp_q[], /* (o) Q15 : Quantized LSP              */
Word16 ana[] /* (o)     : indexes                    */
);

void Get_wegt(const Word16 flsp[], /* (i) Q13 : M LSP parameters  */
Word16 wegt[] /* (o) Q11->norm : M weighting coefficients */
);

void Lsp_encw_reset(Word16 qua_freq_prev[MA_NP][M]);

void Lsp_expand_1(Word16 buf[], /* Q13 */
Word16 gap /* Q13 */
);

void Lsp_expand_2(Word16 buf[], /* Q13 */
Word16 gap /* Q13 */
);

void Lsp_expand_1_2(Word16 buf[], /* Q13 */
Word16 gap /* Q13 */
);

void Lsp_get_quant(const Word16 lspcb1[][M], /* (i) Q13 : first stage LSP codebook      */
const Word16 lspcb2[][M], /* (i) Q13 : Second stage LSP codebook     */
Word16 code0, /* (i)     : selected code of first stage  */
Word16 code1, /* (i)     : selected code of second stage */
Word16 code2, /* (i)     : selected code of second stage */
const Word16 fg[][M], /* (i) Q15 : MA prediction coef.           */
Word16 freq_prev[MA_NP][M], /* (i) Q13 : previous LSP vector           */
Word16 lspq[], /* (o) Q13 : quantized LSP parameters      */
const Word16 fg_sum[] /* (i) Q15 : present MA prediction coef.   */
);

void Lsp_stability(Word16 buf[] /* Q13 */
);

struct lspdec_ctx {
	Word16 freq_prev[MA_NP][M]; /* Q13 */
	/* memory for frame erase operation */
	Word16 prev_ma; /* previous MA prediction coef.*/
	Word16 prev_lsp[M]; /* previous LSP vector         */
};

void D_lsp(struct lspdec_ctx *lspdc, Word16 prm[], /* (i)     : indexes of the selected LSP */
Word16 lsp_q[], /* (o) Q15 : Quantized LSP parameters    */
Word16 erase /* (i)     : frame erase information     */
);

void Lsp_decw_reset(struct lspdec_ctx *lspdc);

void Lsp_prev_compose(const Word16 lsp_ele[M], /* (i) Q13 : LSP vectors                 */
Word16 lsp[M], /* (o) Q13 : quantized LSP parameters    */
const Word16 fg[MA_NP][M], /* (i) Q15 : MA prediction coef.         */
const Word16 freq_prev[MA_NP][M], /* (i) Q13 : previous LSP vector         */
const Word16 fg_sum[M] /* (i) Q15 : present MA prediction coef. */
);

void Lsp_prev_extract(const Word16 lsp[M], /* (i) Q13 : unquantized LSP parameters  */
Word16 lsp_ele[M], /* (o) Q13 : target vector               */
const Word16 fg[MA_NP][M], /* (i) Q15 : MA prediction coef.         */
const Word16 freq_prev[MA_NP][M], /* (i) Q13 : previous LSP vector         */
const Word16 fg_sum_inv[M] /* (i) Q12 : inverse previous LSP vector */
);

void Lsp_prev_update(const Word16 lsp_ele[M], /* (i)   Q13 : LSP vectors           */
Word16 freq_prev[MA_NP][M] /* (i/o) Q13 : previous LSP vectors  */
);

/*-------------------------------*
 * gain VQ constants.            *
 *-------------------------------*/

#define NCODE1_B  3                /* number of Codebook-bit */
#define NCODE2_B  4                /* number of Codebook-bit */
#define NCODE1    (1<<NCODE1_B)    /* Codebook 1 size */
#define NCODE2    (1<<NCODE2_B)    /* Codebook 2 size */
#define NCAN1     4                /* Pre-selecting order for #1 */
#define NCAN2     8                /* Pre-selecting order for #2 */
#define INV_COEF  -17103           /* Q19 */

/*--------------------------------------------------------------------------*
 * gain VQ functions.                                                       *
 *--------------------------------------------------------------------------*/

Word16 Qua_gain(Word16 code[], /* (i) Q13 :Innovative vector.             */
Word16 g_coeff[], /* (i)     :Correlations <xn y1> -2<y1 y1> */
/*            <y2,y2>, -2<xn,y2>, 2<y1,y2> */
Word16 exp_coeff[], /* (i)     :Q-Format g_coeff[]             */
Word16 *gain_pit, /* (o) Q14 :Pitch gain.                    */
Word16 *gain_cod, /* (o) Q1  :Code gain.                     */
Word16 tameflag, /* (i)     : set to 1 if taming is needed  */
Word16 past_qua_en[4]);

void Dec_gain(Word16 index, /* (i)     :Index of quantization.         */
Word16 code[], /* (i) Q13 :Innovative vector.             */
Word16 bfi, /* (i)     :Bad frame indicator            */
Word16 *gain_pit, /* (o) Q14 :Pitch gain.                    */
Word16 *gain_cod, /* (o) Q1  :Code gain.                     */
Word16 past_qua_en[4]);

void Gain_predict(Word16 past_qua_en[], /* (i) Q10 :Past quantized energies        */
Word16 code[], /* (i) Q13 :Innovative vector.             */
Word16 *gcode0, /* (o) Qxx :Predicted codebook gain        */
Word16 *exp_gcode0 /* (o)     :Q-Format(gcode0)               */
);

void Gain_update(Word16 past_qua_en[],/* (i) Q10 :Past quantized energies                  */
Word32 L_gbk12 /* (i) Q13 : gbk1[indice1][1]+gbk2[indice2][1]          */
);

void Gain_update_erasure(Word16 past_qua_en[]/* (i) Q10 :Past quantized energies                   */
);

void Corr_xy2(Word16 xn[], /* (i) Q0  :Target vector.                  */
Word16 y1[], /* (i) Q0  :Adaptive codebook.              */
Word16 y2[], /* (i) Q12 :Filtered innovative vector.     */
Word16 g_coeff[], /* (o) Q[exp]:Correlations between xn,y1,y2 */
Word16 exp_g_coeff[] /* (o)       :Q-format of g_coeff[]         */
);

/*-----------------------*
 * Bitstream function    *
 *-----------------------*/

void prm2bits_ld8k(Word16 prm[], Word16 bits[]);
void bits2prm_ld8k(Word16 bits[], Word16 prm[]);
#define BIT_0     (short)0x007f /* definition of zero-bit in bit-stream      */
#define BIT_1     (short)0x0081 /* definition of one-bit in bit-stream       */
#define SYNC_WORD (short)0x6b21 /* definition of frame erasure flag          */
#define SIZE_WORD (short)80     /* number of speech bits                     */

/*-----------------------------------*
 * Post-filter functions.            *
 *-----------------------------------*/

#define L_H 22     /* size of truncated impulse response of A(z/g1)/A(z/g2) */

#define GAMMAP      16384   /* 0.5               (Q15) */
#define INV_GAMMAP  21845   /* 1/(1+GAMMAP)      (Q15) */
#define GAMMAP_2    10923   /* GAMMAP/(1+GAMMAP) (Q15) */

#define  GAMMA2_PST 18022 /* Formant postfilt factor (numerator)   0.55 Q15 */
#define  GAMMA1_PST 22938 /* Formant postfilt factor (denominator) 0.70 Q15 */

#define  MU       26214   /* Factor for tilt compensation filter   0.8  Q15 */
#define  AGC_FAC  29491   /* Factor for automatic gain control     0.9  Q15 */
#define  AGC_FAC1 (Word16)(32767 - AGC_FAC)    /* 1-AGC_FAC in Q15          */

struct postfilt_ctx {
	/* inverse filtered synthesis (with A(z/GAMMA2_PST))   */
	Word16 res2_buf[PIT_MAX + L_SUBFR];
	Word16 *res2;
	Word16 scal_res2_buf[PIT_MAX + L_SUBFR];
	Word16 *scal_res2;
	/* memory of filter 1/A(z/GAMMA1_PST) */
	Word16 mem_syn_pst[M];

	Word16 mem_pre;

	Word16 past_gain;
};

void Init_Post_Filter(struct postfilt_ctx *postfc);

void Post_Filter(struct postfilt_ctx *postfc, Word16 *syn, /* in/out: synthesis speech (postfiltered is output)    */
Word16 *Az_4, /* input : interpolated LPC parameters in all subframes */
Word16 *T, /* input : decoded pitch lags in all subframes          */
Word16 Vad);

/*--------------------------------------------------------------------------*
 * Constants and prototypes for taming procedure.                           *
 *--------------------------------------------------------------------------*/

#define GPCLIP      15564      /* Maximum pitch gain if taming is needed Q14*/
#define GPCLIP2     481        /* Maximum pitch gain if taming is needed Q9 */
#define GP0999      16383      /* Maximum pitch gain if taming is needed    */
#define L_THRESH_ERR 983040000L /* Error threshold taming 16384. * 60000.   */

void Init_exc_err(Word32 L_exc_err[4]);
void update_exc_err(Word32 L_exc_err[4], Word16 gain_pit, Word16 t0);
Word16 test_err(Word32 L_exc_err[4], Word16 t0, Word16 t0_frac);
/* CNG excitation generation */
void Calc_exc_rand(Word16 cur_gain, /* (i)   :   target sample gain                 */
Word16 *exc, /* (i/o) :   excitation array                   */
Word16 *seed, /* (i)   :   current Vad decision               */
Word32 L_exc_err[4] /* (i)   :   encoder/decoder flag               */
);

/*--------------------------------------------------------------------------*
 * Prototypes for auxiliary functions.                                      *
 *--------------------------------------------------------------------------*/

static inline void Copy(const Word16 x[restrict], /* (i)   : input vector   */
Word16 y[restrict], /* (o)   : output vector  */
Word16 L /* (i)   : vector length  */
) {
    Word16 i;

    for (i = 0; i < L; i++)
        y[i] = x[i];

    return;
};

void Set_zero(Word16 x[], /* (o)    : vector to clear     */
Word16 L /* (i)    : length of vector    */
);

/*-------------------------------*
 * VAD                           *
 *-------------------------------*/
#define     NP            12                  /* Increased LPC order */

struct vad_ctx {
	Word16 MeanLSF[M];
	Word16 Min_buffer[16];
	Word16 Prev_Min, Next_Min, Min;
	Word16 MeanE, MeanSE, MeanSLE, MeanSZC;
	Word16 prev_energy;
	Word16 count_sil, count_update, count_ext;
	Word16 flag, v_flag, less_count;
};

void vad_init(struct vad_ctx *vadc);

void vad(struct vad_ctx *vadc, Word16 rc, Word16 *lsf, Word16 *r_h, Word16 *r_l,
		Word16 exp_R0, Word16 *sigpp, Word16 frm_count, Word16 prev_marker,
		Word16 pprev_marker, Word16 *marker);

/*-------------------------------*
 * DTX                           *
 *-------------------------------*/
#define INIT_SEED       11111
#define A_GAIN0         28672
#define A_GAIN1         4096    /* 32768L - A_GAIN0 */

#define NB_SUMACF       3
#define NB_CURACF       2
#define NB_GAIN         2
#define SIZ_SUMACF      (NB_SUMACF * MP1)
#define SIZ_ACF         (NB_CURACF * MP1)

struct dtx_ctx {
	Word16 lspSid_q[M];
	Word16 pastCoeff[MP1];
	Word16 RCoeff[MP1];
	Word16 sh_RCoeff;
	Word16 Acf[SIZ_ACF];
	Word16 sh_Acf[NB_CURACF];
	Word16 sumAcf[SIZ_SUMACF];
	Word16 sh_sumAcf[NB_SUMACF];
	Word16 ener[NB_GAIN];
	Word16 sh_ener[NB_GAIN];
	Word16 fr_cur;
	Word16 cur_gain;
	Word16 nb_ener;
	Word16 sid_gain;
	Word16 flag_chang;
	Word16 prev_energy;
	Word16 count_fr0;

	/* place levinson coff here because it most use in dtx part */
	Word16 old_A[M + 1];
	Word16 old_rc[2];
};

void Init_Cod_cng(struct dtx_ctx *dtxc);
void Cod_cng(struct dtx_ctx *dtxc, Word16 *exc, /* (i/o) : excitation array                     */
Word16 pastVad, /* (i)   : previous VAD decision                */
Word16 *lsp_old_q, /* (i/o) : previous quantized lsp               */
Word16 *Aq, /* (o)   : set of interpolated LPC coefficients */
Word16 *ana, /* (o)   : coded SID parameters                 */
Word16 freq_prev[MA_NP][M],
/* (i/o) : previous LPS for quantization        */
Word16 *seed, /* (i/o) : random generator seed                */
Word32 L_exc_err[4]);
void Update_cng(struct dtx_ctx *dtxc, Word16 *r_h, /* (i) :   MSB of frame autocorrelation        */
Word16 exp_r, /* (i) :   scaling factor associated           */
Word16 Vad /* (i) :   current Vad decision                */
);

void Levinson(struct dtx_ctx *dtxc, Word16 Rh[], /* (i)     : Rh[m+1] Vector of autocorrelations (msb) */
Word16 Rl[], /* (i)     : Rl[m+1] Vector of autocorrelations (lsb) */
Word16 A[], /* (o) Q12 : A[m]    LPC coefficients  (m = 10)       */
Word16 rc[], /* (o) Q15 : rc[M]   Relection coefficients.          */
Word16 *Err /* (o)     : Residual energy                          */
);

/* SID LSP Quantization */
void Get_freq_prev(Word16 qua_freq_prev[MA_NP][M], Word16 x[MA_NP][M]);
void Update_freq_prev(Word16 qua_freq_prev[MA_NP][M], Word16 x[MA_NP][M]);
void Get_decfreq_prev(struct lspdec_ctx *lspdc, Word16 x[MA_NP][M]);
void Update_decfreq_prev(struct lspdec_ctx *lspdc, Word16 x[MA_NP][M]);

/* Decoder CNG generation */
struct dec_sid_ctx {
	Word16 cur_gain;
	Word16 lspSid[M];
	Word16 sid_gain;
};

void Init_Dec_cng(struct dec_sid_ctx *decsidc);
void Dec_cng(struct dec_sid_ctx *decsidc, Word16 past_ftyp, /* (i)   : past frame type                      */
Word16 sid_sav, /* (i)   : energy to recover SID gain           */
Word16 sh_sid_sav, /* (i)   : corresponding scaling factor         */
Word16 *parm, /* (i)   : coded SID parameters                 */
Word16 *exc, /* (i/o) : excitation array                     */
Word16 *lsp_old, /* (i/o) : previous lsp                         */
Word16 *A_t, /* (o)   : set of interpolated LPC coefficients */
Word16 *seed, /* (i/o) : random generator seed                */
Word16 freq_prev[MA_NP][M]/* (i/o) : previous LPS for quantization        */
);

/*----------------------------------*
 * Main coder and decoder functions *
 *----------------------------------*/
struct cod_ld8a_ctx {
	/* Speech vector */
	Word16 old_speech[L_TOTAL];
	Word16 *speech, *p_window;
	Word16 *new_speech;

	/* Weighted speech vector */
	Word16 old_wsp[L_FRAME + PIT_MAX];
	Word16 *wsp;

	/* Excitation vector */
	Word16 old_exc[L_FRAME + PIT_MAX + L_INTERPOL];
	Word16 *exc;

	/* Lsp (Line spectral pairs) */
	Word16 lsp_old[M];
	Word16 lsp_old_q[M];

	/* Filter's memory */
	Word16 mem_w0[M], mem_w[M], mem_zero[M];
	Word16 sharp;

	/* Gain predictor */
	Word16 past_qua_en[4];

	/* For G.729B */
	/* DTX variables */
	Word16 pastVad;
	Word16 ppastVad;
	Word16 seed;

	struct vad_ctx vadc;

	struct prepost_process_ctx prepc;

	struct dtx_ctx dtxc;

	Word16 qua_freq_prev[MA_NP][M];
	/* taming */
	Word32 L_exc_err[4];
};

void Init_Coder_ld8a(struct cod_ld8a_ctx *codc);

void Coder_ld8a(struct cod_ld8a_ctx *codc, Word16 ana[], /* output  : Analysis parameters */
Word16 frame, /* input   : frame counter       */
Word16 vad_enable /* input   : VAD enable flag     */
);

struct dec_ld8a_ctx {
	/* Excitation vector */
	Word16 old_exc[L_FRAME + PIT_MAX + L_INTERPOL];
	Word16 *exc;

	/* Gain predictor */
	Word16 past_qua_en[4];

	/* Lsp (Line spectral pairs) */
	Word16 lsp_old[M];

	/* Filter's memory */
	Word16 mem_syn[M];

	Word16 sharp; /* pitch sharpening of previous frame */
	Word16 old_T0; /* integer delay of previous frame    */
	Word16 gain_code; /* Code gain                          */
	Word16 gain_pitch; /* Pitch gain                         */

	/* for G.729B */
	Word16 seed_fer;
	/* CNG variables */
	Word16 past_ftyp;
	Word16 seed;
	Word16 sid_sav, sh_sid_sav;

	struct lspdec_ctx lspdc;

	struct dec_sid_ctx decsidc;

	struct postfilt_ctx postfc;

	struct prepost_process_ctx postprc;
};

void Init_Decod_ld8a(struct dec_ld8a_ctx *decc);

void Decod_ld8a(struct dec_ld8a_ctx *decc, Word16 parm[], /* (i)   : vector of synthesis parameters
 parm[0] = bad frame indicator (bfi)  */
Word16 synth[], /* (o)   : synthesis speech                     */
Word16 A_t[], /* (o)   : decoded LP filter in 2 subframes     */
Word16 *T2, /* (o)   : decoded pitch lag in 2 subframes     */
Word16 *Vad /* (o)   : frame type                           */
);
