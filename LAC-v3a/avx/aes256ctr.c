/*
  crypto_stream_aes256ctr
  based heavily on public-domain code by Romain Dolbeau
  Different handling of nonce+counter than original version
  using separated 96-bit nonce and internal 32-bit counter, starting from zero
  Public Domain
*/

#include <stdint.h>
#include <emmintrin.h>
#include <immintrin.h>
#include "aes256ctr.h"

static inline void aesni_encrypt8(unsigned char *out,
                                  __m128i *n,
                                  const __m128i rkeys[16])
{
  __m128i nv0;
  __m128i nv1;
  __m128i nv2;
  __m128i nv3;
  __m128i nv4;
  __m128i nv5;
  __m128i nv6;
  __m128i nv7;

  /* Load current counter value */
  __m128i nv0i = _mm_load_si128(n);

  /* Increase counter in 8 consecutive blocks */
  nv0 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(0,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv1 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(1,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv2 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(2,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv3 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(3,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv4 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(4,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv5 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(5,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv6 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(6,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));
  nv7 = _mm_shuffle_epi8(_mm_add_epi32(nv0i, _mm_set_epi64x(7,0)), _mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));

  /* Write counter for next iteration, increased by 8 */
  _mm_store_si128(n, _mm_add_epi32(nv0i, _mm_set_epi64x(8,0)));

  /* Actual AES encryption, 8x interleaved */
  __m128i temp0 = _mm_xor_si128(nv0, rkeys[0]);
  __m128i temp1 = _mm_xor_si128(nv1, rkeys[0]);
  __m128i temp2 = _mm_xor_si128(nv2, rkeys[0]);
  __m128i temp3 = _mm_xor_si128(nv3, rkeys[0]);
  __m128i temp4 = _mm_xor_si128(nv4, rkeys[0]);
  __m128i temp5 = _mm_xor_si128(nv5, rkeys[0]);
  __m128i temp6 = _mm_xor_si128(nv6, rkeys[0]);
  __m128i temp7 = _mm_xor_si128(nv7, rkeys[0]);

  for (int i = 1; i < 14; i++) {
    temp0 =  _mm_aesenc_si128(temp0, rkeys[i]);
    temp1 =  _mm_aesenc_si128(temp1, rkeys[i]);
    temp2 =  _mm_aesenc_si128(temp2, rkeys[i]);
    temp3 =  _mm_aesenc_si128(temp3, rkeys[i]);
    temp4 =  _mm_aesenc_si128(temp4, rkeys[i]);
    temp5 =  _mm_aesenc_si128(temp5, rkeys[i]);
    temp6 =  _mm_aesenc_si128(temp6, rkeys[i]);
    temp7 =  _mm_aesenc_si128(temp7, rkeys[i]);
  }

  temp0 = _mm_aesenclast_si128(temp0, rkeys[14]);
  temp1 = _mm_aesenclast_si128(temp1, rkeys[14]);
  temp2 = _mm_aesenclast_si128(temp2, rkeys[14]);
  temp3 = _mm_aesenclast_si128(temp3, rkeys[14]);
  temp4 = _mm_aesenclast_si128(temp4, rkeys[14]);
  temp5 = _mm_aesenclast_si128(temp5, rkeys[14]);
  temp6 = _mm_aesenclast_si128(temp6, rkeys[14]);
  temp7 = _mm_aesenclast_si128(temp7, rkeys[14]);

  /* Write results */
  _mm_storeu_si128((__m128i*)(out+  0), temp0);
  _mm_storeu_si128((__m128i*)(out+ 16), temp1);
  _mm_storeu_si128((__m128i*)(out+ 32), temp2);
  _mm_storeu_si128((__m128i*)(out+ 48), temp3);
  _mm_storeu_si128((__m128i*)(out+ 64), temp4);
  _mm_storeu_si128((__m128i*)(out+ 80), temp5);
  _mm_storeu_si128((__m128i*)(out+ 96), temp6);
  _mm_storeu_si128((__m128i*)(out+112), temp7);
}

void aes256ctr_init(aes256ctr_ctx *state,
                    const unsigned char *key,
                    uint16_t nonce)
{
  __m128i key0 = _mm_loadu_si128((__m128i *)(key+0));
  __m128i key1 = _mm_loadu_si128((__m128i *)(key+16));
  __m128i temp0, temp1, temp2, temp4;
  int idx = 0;

  state->n = _mm_set_epi64x(0, (uint64_t)nonce << 48);

  state->rkeys[idx++] = key0;
  temp0 = key0;
  temp2 = key1;
  temp4 = _mm_setzero_si128();

#define BLOCK1(IMM)                                                     \
  temp1 = _mm_aeskeygenassist_si128(temp2, IMM);                        \
  state->rkeys[idx++] = temp2;                                          \
  temp4 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(temp4), _mm_castsi128_ps(temp0), 0x10));  \
  temp0 = _mm_xor_si128(temp0, temp4);                                  \
  temp4 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(temp4), _mm_castsi128_ps(temp0), 0x8c));  \
  temp0 = _mm_xor_si128(temp0, temp4);                                  \
  temp1 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(temp1), _mm_castsi128_ps(temp1), 0xff));  \
  temp0 = _mm_xor_si128(temp0, temp1)

#define BLOCK2(IMM)                                                     \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                        \
  state->rkeys[idx++] = temp0;                                          \
  temp4 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(temp4), _mm_castsi128_ps(temp2), 0x10));  \
  temp2 = _mm_xor_si128(temp2, temp4);                                  \
  temp4 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(temp4), _mm_castsi128_ps(temp2), 0x8c));  \
  temp2 = _mm_xor_si128(temp2, temp4);                                  \
  temp1 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(temp1), _mm_castsi128_ps(temp1), 0xaa));  \
  temp2 = _mm_xor_si128(temp2, temp1)

  BLOCK1(0x01);
  BLOCK2(0x01);

  BLOCK1(0x02);
  BLOCK2(0x02);

  BLOCK1(0x04);
  BLOCK2(0x04);

  BLOCK1(0x08);
  BLOCK2(0x08);

  BLOCK1(0x10);
  BLOCK2(0x10);

  BLOCK1(0x20);
  BLOCK2(0x20);

  BLOCK1(0x40);
  state->rkeys[idx++] = temp0;
}

void aes256ctr_select(aes256ctr_ctx *state, uint16_t nonce) {
  state->n = _mm_set_epi64x(0, (uint64_t)nonce << 48);
}

void aes256ctr_squeezeblocks(unsigned char *out,
                             unsigned long long nblocks,
                             aes256ctr_ctx *state)
{
  unsigned long long i;

  for(i=0;i<nblocks;i++) {
    aesni_encrypt8(out, &state->n, state->rkeys);
    out += 128;
  }
}

void aes256ctr_prf(unsigned char *out,
                   unsigned long long outlen,
                   const unsigned char *seed,
                   unsigned char nonce)
{
  unsigned int i;
  unsigned char buf[128];
  aes256ctr_ctx state;

  aes256ctr_init(&state, seed, (uint16_t)nonce << 8);

  while(outlen >= 128) {
    aesni_encrypt8(out, &state.n, state.rkeys);
	out+=128;
    outlen -= 128;
  }

  if(outlen) {
    aesni_encrypt8(buf, &state.n, state.rkeys);
    for(i=0;i<outlen;i++)
      out[i] = buf[i];
  }
}
