#include <g729itu.h>
#include <string.h>
#include "../Sources/typedef.h"
#include "../Sources/ld8a.h"
#include "../Sources/basic_op.h"

#include "tstseq1_bin.h"
#include "tstseq1a_bit.h"
#include "tstseq1a_out.h"
#include "tstseq2_bin.h"
#include "tstseq2a_bit.h"
#include "tstseq2a_out.h"
#include "tstseq3_bin.h"
#include "tstseq3a_bit.h"
#include "tstseq3a_out.h"
#include "tstseq4_bin.h"
#include "tstseq4a_bit.h"
#include "tstseq4a_out.h"

#ifdef VECTORS_TEST

void g729itu_prm_decode(struct g729itu_vocoder *voc, Word16 prm[PRM_SIZE + 1],
		int16_t signal[]);
void g729itu_prm_encode(struct g729itu_vocoder *voc, int16_t inputFrame[],
		Word16 parm[PRM_SIZE + 1]);
int read_frame(const uint8_t *bitvector, Word16 *parm);

static int coder(struct g729itu_vocoder *voc, const uint8_t *input,
		uint32_t input_size, uint8_t *output) {
	int read_ind = 0;
	int write_ind = 0;
	int ser_size;
	Word16 prm[PRM_SIZE + 1]; /* Analysis parameters + frame type      */
	Word16 serial[SERIAL_SIZE]; /* Output bitstream buffer               */

	while (input_size >= L_FRAME * sizeof(Word16)) {
		g729itu_prm_encode(voc, (int16_t *) (input + read_ind), prm);
		read_ind += L_FRAME * sizeof(Word16);
		input_size -= L_FRAME * sizeof(Word16);

		prm2bits_ld8k(prm, serial);
		ser_size = (serial[1] + (Word16) 2) * sizeof(Word16);
		memcpy(output + write_ind, serial, ser_size);
		write_ind += ser_size;
	}

	return write_ind;
}

static int decoder(struct g729itu_vocoder *voc, const uint8_t *input,
		uint32_t input_size, uint8_t *output) {
	int read_ind = 0;
	int write_ind = 0;
	Word16 prm[PRM_SIZE + 2]; /* Analysis parameters + frame type      */
	Word16 serial[L_FRAME]; /* Output bitstream buffer               */

	while (read_ind < input_size) {
		read_ind += read_frame(input + read_ind, prm);

		g729itu_prm_decode(voc, prm, serial);

		memcpy(output + write_ind, serial, sizeof(Word16) * L_FRAME);
		write_ind += sizeof(Word16) * L_FRAME;
	}

	return write_ind;
}

static uint8_t perform_test(struct g729itu_vocoder *voc, uint8_t *temp_buffer,
		const unsigned char *bin, size_t bin_size, const unsigned char *bit,
		size_t bit_size, const unsigned char *out, size_t out_size) {
	uint8_t test_ret = 0;
	int ret;
	g729itu_init(voc, 1);
	Overflow = Carry = 0;
	ret = coder(voc, bin, bin_size, temp_buffer);
	if (ret != bit_size) {
		test_ret |= 1;
	} else {
		ret = memcmp(temp_buffer, bit, bit_size);
		if (ret) {
			test_ret |= 2;
		}
	}

	g729itu_init(voc, 1);
	ret = decoder(voc, bit, bit_size, temp_buffer);
	if (ret != out_size) {
		test_ret |= 4;
	} else {
		ret = memcmp(temp_buffer, out, out_size);
		if (ret) {
			test_ret |= 8;
		}
	}

	return test_ret;
}

size_t g729_test_get_test_size(void) {
	return tstseq1_bin_size;
}

uint16_t g729_test(struct g729itu_vocoder *voc, uint8_t *test_buffer) {
	uint16_t res;

	res = perform_test(voc, test_buffer, tstseq1_bin, tstseq1_bin_size,
			tstseq1a_bit, tstseq1a_bit_size, tstseq1a_out, tstseq1a_out_size);
	res <<= 4;

	res |= perform_test(voc, test_buffer, tstseq2_bin, tstseq2_bin_size,
			tstseq2a_bit, tstseq2a_bit_size, tstseq2a_out, tstseq2a_out_size);
	res <<= 4;

	res |= perform_test(voc, test_buffer, tstseq3_bin, tstseq3_bin_size,
			tstseq3a_bit, tstseq3a_bit_size, tstseq3a_out, tstseq3a_out_size);
	res <<= 4;

	res |= perform_test(voc, test_buffer, tstseq4_bin, tstseq4_bin_size,
			tstseq4a_bit, tstseq4a_bit_size, tstseq4a_out, tstseq4a_out_size);

	return res;
}
#endif
