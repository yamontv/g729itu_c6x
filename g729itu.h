#ifndef HEADERS_G729ITU_H_
#define HEADERS_G729ITU_H_

#include <stdlib.h>
#include <stdint.h>

/* comment this to not include test data (big) */
//#define VECTORS_TEST

/* front declaration */
struct g729itu_vocoder;

/* get size of vocoder in bytes
 * you have to memalign space by yourself
 * ForEx:
 * struct g729itu_vocoder *voc = memalign(128, g729itu_get_vocoder_size());
 * */
size_t g729itu_get_vocoder_size(void);

/* initialize vocoder context
 * run it once before using encode decode
 * enableVAD - 0 means don't use AnexB
 */
void g729itu_init(struct g729itu_vocoder *voc, uint8_t enableVAD);

/* encode 10ms of linear audio
 * inputFrame - input 80 linear samples (10ms)
 * bitStream - output encoded data ready to be transmitted over RTP
 * bitStreamLength - output length of bitStream data in bytes
 * ForEx:
 * g729itu_encode(voc, input_linear, rtp_payload_out, &rtp_payload_len);
 */
void g729itu_encode(struct g729itu_vocoder *ctx, const int16_t inputFrame[],
		uint8_t bitStream[], uint8_t *bitStreamLength);

/* decode 10ms of encoded data
 * bitStream - input encoded data from RTP stream
 * bitStreamLength - input encoded data length
 * signal - output decoded 80 linear samples (10ms)
 * ForEx:
 * g729itu_decode(voc, rtp_payload_in, rtp_payload_len, output_linear);
 */
void g729itu_decode(struct g729itu_vocoder *ctx, const uint8_t bitStream[],
		uint8_t bitStreamLength, int16_t signal[]);

/* Debug output from codec
 * MUST be implement!
 * you can create empty function if you want
 * msg - null-terminated debug string from codec*/
void g729_itu_debug_out(const char *msg);

#ifdef VECTORS_TEST
/* get size of test_buffer in bytes
 * you have to alloc space by yourself
 * no align requirements
 * ForEx:
 * uint8_t *test_buffer = malloc(g729_test_get_test_size());
 * */
size_t g729_test_get_test_size(void);

/* perform vector test
 * return 0 if success
 * ForEx:
 * struct g729itu_vocoder *voc = memalign(128, g729itu_get_vocoder_size());
 * uint8_t *test_buffer = malloc(g729_test_get_test_size());
 * ret = g729_test(voc, test_buffer);
 */
uint16_t g729_test(struct g729itu_vocoder *voc, uint8_t *test_buffer);
#endif

#endif /* HEADERS_G729ITU_H_ */
