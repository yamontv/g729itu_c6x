/* ITU-T G.729 Software Package Release 2 (November 2006) */
/*
 ITU-T G.729A Annex B     ANSI-C Source Code
 Version 1.3    Last modified: August 1997
 Copyright (c) 1996, France Telecom, Rockwell International,
 Universite de Sherbrooke.
 All rights reserved.
 */

#define         TRUE 1
#define         FALSE 0
#define         sqr(a)  ((a)*(a))
#define         R_LSFQ 10

void lsfq_noise(Word16 *lsp, Word16 *lspq, Word16 freq_prev[MA_NP][M],
		Word16 *ana);
void sid_lsfq_decode(Word16 *index, /* (i) : quantized indices    */
Word16 *lspq, /* (o) : quantized lsp vector */
Word16 freq_prev[MA_NP][M] /* (i) : memory of predictor  */
);

void Qua_Sidgain(Word16 *ener, /* (i)   array of energies                   */
Word16 *sh_ener, /* (i)   corresponding scaling factors       */
Word16 nb_ener, /* (i)   number of energies or               */
Word16 *enerq, /* (o)   decoded energies in dB              */
Word16 *idx /* (o)   SID gain quantization index         */
);

