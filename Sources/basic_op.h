#include <c6x.h>
#include <gsm.h>

#define g_round round

#define get_sat() (_extu(CSR, 22, 31))

#define clr_sat() (CSR = _clr(CSR, 9, 9))

static inline Word16 Random(Word16 *seed) {

    /* seed = seed*31821 + 13849; */
    *seed = extract_l(L_add(L_shr(L_mult(*seed, 31821), 1), 13849L));

    return (*seed);
};
