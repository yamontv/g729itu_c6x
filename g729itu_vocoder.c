#include <g729itu.h>
#include <string.h>
#include "Sources/typedef.h"
#include "Sources/ld8a.h"
#include "Sources/tab_ld8a.h"
#include "Sources/sid.h"

struct g729itu_vocoder {
	struct cod_ld8a_ctx coder;
	struct dec_ld8a_ctx decoder;

	/* coder data */
	Word16 cprm[PRM_SIZE + 1];
	Word32 cframe;

	/* decoder data */
	Word16 dsynth_buf[L_FRAME + M], *synth; /* Synthesis                   */
	Word16 dparm[PRM_SIZE + 2]; /* Synthesis parameters        */
	Word16 dAz_dec[MP1 * 2]; /* Decoded Az for post-filter  */
	Word16 dT2[2]; /* Pitch lag for 2 subframes   */

	int vad_enable;
};

size_t g729itu_get_vocoder_size(void) {
	return sizeof(struct g729itu_vocoder);
}

void g729itu_init(struct g729itu_vocoder *voc, uint8_t enableVAD) {
	memset(voc, 0, sizeof(struct g729itu_vocoder));
	voc->vad_enable = enableVAD;

	/* init coder */
	Init_Pre_Process(&voc->coder.prepc);
	Init_Coder_ld8a(&voc->coder);
	Init_Cod_cng(&voc->coder.dtxc);

	/* init decoder */
	Init_Decod_ld8a(&voc->decoder);
	Init_Post_Filter(&voc->decoder.postfc);
	Init_Post_Process(&voc->decoder.postprc);
	voc->synth = voc->dsynth_buf + M;

	/* for G.729b */
	Init_Dec_cng(&voc->decoder.decsidc);
}

static const unsigned short bit_mask[16] = { 0x0001, 0x0003, 0x0007, 0x000F,
		0x001F, 0x003F, 0x007F, 0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF,
		0x3FFF, 0x7FFF, 0xFFFF };

static inline void g729_prm2frm(uint8_t *frame, uint8_t *frm_len,
		Word16 param[]) {
	Word16 bit_num = 0;
	Word16 i, tmp_word = 0;

	switch (param[0]) {
	case 0:
	default:
		/* not transmitted */
		*frm_len = 0;
		break;

	case 1:
		for (i = 0; i < PRM_SIZE; i++) {
			tmp_word = (tmp_word << bitsno[i]) | param[i + 1];
			bit_num += bitsno[i];
			while (bit_num >= 8) {
				bit_num = bit_num - 8;
				*frame++ = (tmp_word >> bit_num) & 0xFF;
				if (bit_num > 0)
					tmp_word &= bit_mask[bit_num - 1];
				else
					tmp_word = 0;
			}
		}
		*frm_len = 10;
		break;

	case 2:
		for (i = 0; i < 4; i++) {
			tmp_word = (tmp_word << bitsno2[i]) | param[i + 1];
			bit_num += bitsno2[i];
			while (bit_num >= 8) {
				bit_num = bit_num - 8;
				*frame++ = (tmp_word >> bit_num) & 0xFF;
				if (bit_num > 0)
					tmp_word &= bit_mask[bit_num - 1];
				else
					tmp_word = 0;
			}
		}
		if ((bit_num != 0) && (bit_num < 8))
			*frame++ = tmp_word << (8 - bit_num);
		/* SID frame */
		*frm_len = 2;
		break;
	}
	return;
}

static inline void g729_frm2prm(const uint8_t *frame, uint8_t frm_len,
		Word16 param[]) {
	Word16 bit_num = 0;
	Word16 i, tmp_word = 0;

	/* Decode the frame. */
	if (frm_len == 10) {
		/* 80 bits = 10 bytes (speech data) */
		param[1] = 1;

		for (i = 0; i < PRM_SIZE; i++) {
			while (bit_num < bitsno[i]) {
				tmp_word = (tmp_word << 8) | *frame++;
				bit_num = bit_num + 8;
			}
			bit_num = bit_num - bitsno[i];
			param[i + 2] = tmp_word >> bit_num;
			if (bit_num > 0)
				tmp_word &= bit_mask[bit_num - 1];
			else
				tmp_word = 0;
		}

		/* Check parity and put 1 in param[5] if parity error */
		param[5] = Check_Parity_Pitch(param[4], param[5]);
	} else if (frm_len == 2) {
		/* 15 bits = 2 bytes (SID data) */
		param[1] = 2;

		/* TODO: SS (VAD) doesn't work smoothly!
		 * After 1 minute it stalls somewhere. */
		for (i = 0; i < 4; i++) {
			while (bit_num < bitsno2[i]) {
				tmp_word = (tmp_word << 8) | *frame++;
				bit_num = bit_num + 8;
			}
			bit_num = bit_num - bitsno2[i];
			param[i + 2] = tmp_word >> bit_num;
			if (bit_num > 0)
				tmp_word &= bit_mask[bit_num - 1];
			else
				tmp_word = 0;
		}
	} else {
		param[1] = 0;
	}

	param[0] = 0; /* No frame erasure */

	return;
}

void g729itu_encode(struct g729itu_vocoder *voc, const int16_t inputFrame[],
		uint8_t bitStream[], uint8_t *bitStreamLength) {
	if (voc->cframe == 32767)
		voc->cframe = 256;
	else
		voc->cframe++;

	Copy(inputFrame, voc->coder.new_speech, L_FRAME);

	Pre_Process(&voc->coder.prepc, voc->coder.new_speech);
	Coder_ld8a(&voc->coder, voc->cprm, voc->cframe, voc->vad_enable);
	g729_prm2frm(bitStream, bitStreamLength, voc->cprm);
}

void g729itu_decode(struct g729itu_vocoder *voc, const uint8_t bitStream[],
		uint8_t bitStreamLength, int16_t signal[]) {
	Word16 Vad;

	g729_frm2prm(bitStream, bitStreamLength, voc->dparm);
	Decod_ld8a(&voc->decoder, voc->dparm, voc->synth, voc->dAz_dec, voc->dT2,
			&Vad);
	Post_Filter(&voc->decoder.postfc, voc->synth, voc->dAz_dec, voc->dT2, Vad); /* Post-filter */
	Post_Process(&voc->decoder.postprc, voc->synth);

	Copy(voc->synth, signal, L_FRAME);
}

#ifdef VECTORS_TEST
/* only for itu test */
void g729itu_prm_decode(struct g729itu_vocoder *voc, Word16 prm[PRM_SIZE + 1],
		int16_t signal[]) {
	Word16 Vad;

	Decod_ld8a(&voc->decoder, prm, voc->synth, voc->dAz_dec, voc->dT2, &Vad);
	Post_Filter(&voc->decoder.postfc, voc->synth, voc->dAz_dec, voc->dT2, Vad); /* Post-filter */
	Post_Process(&voc->decoder.postprc, voc->synth);

	Copy(voc->synth, signal, L_FRAME);
}

void g729itu_prm_encode(struct g729itu_vocoder *voc, int16_t inputFrame[],
		Word16 parm[PRM_SIZE + 1]) {
	if (voc->cframe == 32767)
		voc->cframe = 256;
	else
		voc->cframe++;

	Copy(inputFrame, voc->coder.new_speech, L_FRAME);

	Pre_Process(&voc->coder.prepc, voc->coder.new_speech);
	Coder_ld8a(&voc->coder, parm, voc->cframe, voc->vad_enable);
}
#endif
