#ifndef MY_BWT_H_
#define MY_BWT_H_

#include <stdint.h>

typedef uint64_t bwtint_t;
typedef uint8_t ubyte_t;

typedef struct {
	bwtint_t seq_len, bwt_size, n_occ;
	bwtint_t primary;
	uint32_t *bwt;
	bwtint_t *sa, L2[5];
	uint32_t *occ;
	uint32_t cnt_table[256];
} bwtl_t;

#define bwtl_B0(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

#ifdef __cplusplus
extern "C" {
#endif

	int my_bwt_match_exact(const bwtl_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end);
	bwtl_t *my_restore_bwt(const char *fn);
	//bwtl_t *bwtl_seq2bwtl(int len, const uint8_t *seq);
	bwtl_t *bwtl_seq2bwtl();
	inline uint32_t bwtl_occ(const bwtl_t *bwt, uint32_t k, uint8_t c);
	inline void bwtl_occ4(const bwtl_t *bwt, uint32_t k, uint32_t cnt[4]);
	inline void bwtl_2occ4(const bwtl_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4]);
	void bwtl_destroy(bwtl_t *bwt);

#ifdef __cplusplus
}
#endif

#endif
