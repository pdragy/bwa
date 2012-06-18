#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "my_bwt.h"
//#include "bwt_lite.h"
#include "utils.h"
#include "utils.c"

int is_sa(const uint8_t *T, uint32_t *SA, int n);
int is_bwt(uint8_t *T, int n);

int main() {

	const uint8_t *sequence = (uint8_t*)calloc(sizeof(uint8_t*),6);
	//sequence[0] = "A";
	int len = 6;
	//bwtl_t * my_bwt = bwtl_seq2bwtl(len, sequence);
	bwtl_t * my_bwt = bwtl_seq2bwtl();
	
	fprintf(stderr, "restoring bwt: bwt_size=%u, bwt primary=%u, bwt seq_len=%u, L2[0]=%u, L2[1]=%u, L2[2]=%u, L2[3]=%u, L2[4]=%u\n",my_bwt->bwt_size, my_bwt->primary, my_bwt->seq_len,my_bwt->L2[0],my_bwt->L2[1],my_bwt->L2[2],my_bwt->L2[3],my_bwt->L2[4]);
	return 0;
}

int my_bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int my_bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}



bwtl_t *my_restore_bwt(const char *fn)
{
	fprintf(stderr, "restoring bwt\n");
	bwtl_t *b;
	FILE *fp;
	b = (bwtl_t*)calloc(1, sizeof(bwtl_t));
	//fp = fopen(fn, "rb");
	fp = xopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	b->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	b->bwt = (uint32_t*)calloc(b->bwt_size, 4);
	fseek(fp, 0, SEEK_SET);
	fread(&b->primary, sizeof(bwtint_t), 1, fp);
	fread(b->L2+1, sizeof(bwtint_t), 4, fp);
	fread(b->bwt, 4, b->bwt_size, fp);
	b->seq_len = b->L2[4];
	fclose(fp);
	fprintf(stderr, "restoring bwt: bwt_size=%u, bwt primary=%u, bwt seq_len=%u, L2[0]=%u, L2[1]=%u, L2[2]=%u, L2[3]=%u, L2[4]=%u\n",b->bwt_size, b->primary, b->seq_len,b->L2[0],b->L2[1],b->L2[2],b->L2[3],b->L2[4]);
	//bwt_gen_cnt_table(bwt);
	fprintf(stderr, "returning bwt\n");
	return b;
}


//bwtl_t *bwtl_seq2bwtl(int len, const uint8_t *seq)
bwtl_t *bwtl_seq2bwtl()
{
	bwtl_t *b;
	int i;
	//b = (bwtl_t*)calloc(1, sizeof(bwtl_t));
	//b->seq_len = len;

	//b->bwt_size = (len + 15) / 16;
	//b->bwt = (uint32_t*)calloc(b->bwt_size, 4);
	//b->bwt = my_restore_bwt("data/reference/test2.fa.bwt");
	b = my_restore_bwt("data/reference/test2.fa.bwt");
	fprintf(stderr, "restoring bwt: bwt_size=%u, bwt primary=%u, bwt seq_len=%u, L2[0]=%u, L2[1]=%u, L2[2]=%u, L2[3]=%u, L2[4]=%u\n",b->bwt_size, b->primary, b->seq_len,b->L2[0],b->L2[1],b->L2[2],b->L2[3],b->L2[4]);

	int len = b->seq_len;
	//don't need this because we are calculating bwt from file
	/*
	{ // calculate b->bwt
		uint8_t *s;
		b->sa = (uint32_t*)calloc(len + 1, 4);
		is_sa(seq, b->sa, len);
		s = (uint8_t*)calloc(len + 1, 1);
		for (i = 0; i <= len; ++i) {
			if (b->sa[i] == 0) b->primary = i;
			else s[i] = seq[b->sa[i] - 1];
		}
		for (i = b->primary; i < len; ++i) s[i] = s[i + 1];
		b->bwt_size = (len + 15) / 16;
		b->bwt = (uint32_t*)calloc(b->bwt_size, 4);
		for (i = 0; i < len; ++i)
			b->bwt[i>>4] |= s[i] << ((15 - (i&15)) << 1);
		free(s);
	}
	*/

	{ // calculate b->occ
		uint32_t c[4];
		b->n_occ = (len + 15) / 16 * 4;
		b->occ = (uint32_t*)calloc(b->n_occ, 4);
		memset(c, 0, 16);
		for (i = 0; i < len; ++i) {
			if (i % 16 == 0)
				memcpy(b->occ + (i/16) * 4, c, 16);
			++c[bwtl_B0(b, i)];
		}
		memcpy(b->L2+1, c, 16);
		for (i = 2; i < 5; ++i) b->L2[i] += b->L2[i-1];
	}

	{ // generate cnt_table
		for (i = 0; i != 256; ++i) {
			u_int32_t j, x = 0;
			for (j = 0; j != 4; ++j)
				x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
			b->cnt_table[i] = x;
		}
	}
	return b;
}
inline uint32_t bwtl_occ(const bwtl_t *bwt, uint32_t k, uint8_t c)
{
	uint32_t n, b;
	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (uint32_t)(-1)) return 0;
	if (k >= bwt->primary) --k; // because $ is not in bwt
	n = bwt->occ[k/16<<2|c];
	b = bwt->bwt[k/16] & ~((1U<<((15-(k&15))<<1)) - 1);
	n += (bwt->cnt_table[b&0xff] + bwt->cnt_table[b>>8&0xff]
		  + bwt->cnt_table[b>>16&0xff] + bwt->cnt_table[b>>24]) >> (c<<3) & 0xff;
	if (c == 0) n -= 15 - (k&15); // corrected for the masked bits
	return n;
}
inline void bwtl_occ4(const bwtl_t *bwt, uint32_t k, uint32_t cnt[4])
{
	uint32_t x, b;
	if (k == (uint32_t)(-1)) {
		memset(cnt, 0, 16);
		return;
	}
	if (k >= bwt->primary) --k; // because $ is not in bwt
	memcpy(cnt, bwt->occ + (k>>4<<2), 16);
	b = bwt->bwt[k>>4] & ~((1U<<((~k&15)<<1)) - 1);
	x = bwt->cnt_table[b&0xff] + bwt->cnt_table[b>>8&0xff]
		+ bwt->cnt_table[b>>16&0xff] + bwt->cnt_table[b>>24];
	x -= 15 - (k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}
inline void bwtl_2occ4(const bwtl_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
	bwtl_occ4(bwt, k, cntk);
	bwtl_occ4(bwt, l, cntl);
}
void bwtl_destroy(bwtl_t *bwt)
{
	if (bwt) {
		free(bwt->occ); free(bwt->bwt); free(bwt->sa);
		free(bwt);
	}
}
