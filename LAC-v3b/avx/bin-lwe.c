#include <immintrin.h>
#include <string.h>
#include "bin-lwe.h"
#include "rand.h"
#include "lac_param.h"

#ifdef LINUX
#define ALIGNED(x) __attribute__((aligned(x)))
#endif

#ifdef WIN
#define ALIGNED(x) __declspec(align(x))
#endif

//generate the public parameter a from seed
int gen_a(unsigned char *a,  const unsigned char *seed)
{
	//check the pointers
	if(a==NULL || seed==NULL)
	{
		return 1;
	}
	
	pseudo_random_bytes(a,DIM_N,seed);
	//use first 32 bytes of a as seed to generate buf data for replace
	
	return 0;
}
 
//generate the small random vector for secret and error, with fixed hamming weight
//use for e,e1,e2
int gen_e(unsigned char *e,  unsigned char *seed)
{
	
	if(e==NULL)
	{
		return 1;
	}
	
	int i;
	uint16_t buf[NUM_ONE*2];
	gen_r((unsigned char *)buf,seed);
	memset(e,0,DIM_N);
	for(i=0;i<NUM_ONE;i++)
	{
		e[buf[i]]=1;
		e[buf[i+NUM_ONE]]=Q-1;
	}
	
	return 0;
}
 
//for r,s
int gen_r(unsigned char *r,  unsigned char *seed)
{
	
	if(r==NULL)
	{
		return 1;
	}
	
	int i,p;
	uint16_t tmp;
	uint16_t  r_buf[DIM_N],index[SAMPLE_LEN],tmp_index,index_mk;
	uint16_t mk=DIM_N-1;
	unsigned int mask_p,loop=SAMPLE_LEN;

	//init r to be 1,2,3,4,5
	for(i=0;i<DIM_N;i++)
	{
		r_buf[i]=i;
	}
	
	
	p=0;
	while(p<NUM_ONE*2)
	{
		pseudo_random_bytes((unsigned char *)index,SAMPLE_LEN*2,seed);
	    //shuffle
		for(i=0;i<loop;i++)
		{
			//check index
			index_mk=index[i]&mk;
			tmp_index=index_mk>=p ? index_mk: p;
			mask_p=index_mk>=p ? 1:0;
		    //swap
			tmp=r_buf[tmp_index];
			r_buf[tmp_index]=r_buf[p];
			r_buf[p]=tmp;
			//move to the next
			p+=mask_p;
		}
		//update seed
		if(p<NUM_ONE*2)
		{
			memcpy(seed,(unsigned char *)index,SEED_LEN);
		}
	}
	//compy the first NUM_ONE positions to r
	memcpy(r,r_buf,NUM_ONE*2*sizeof(uint16_t));
	
	return 0;
}



// poly_mul  b=as
int poly_mul(const unsigned char *a, const unsigned char *s, unsigned char *b, unsigned int vec_num)
{
	int i,j;
	int8_t  ALIGNED(32) v[DIM_N*2];
	uint16_t *s_one,*s_minusone;
	int8_t *v1_p,*v2_p;
	
	int8_t  ALIGNED(32) sum1[DIM_N],sum2[DIM_N];
	__m256i tmp_sum1,  tmp_v1;
	__m256i tmp_sum2,  tmp_v2;
	
	//extend a to -a_0 -a_1, ...-a_n-1,a_0,a_1,...a_n-1
	for(i=0;i<DIM_N;i++)
	{
		v[i]=-a[i];
	}	
	memcpy(v+DIM_N,a,DIM_N);
	//construt s

	memset(sum1,0,DIM_N);
	memset(sum2,0,DIM_N);
	s_one=(uint16_t*)s;
	s_minusone=(uint16_t*)s+NUM_ONE;
	
	for(i=0;i<NUM_ONE;i++)
	{
		v1_p=(v+DIM_N-s_one[i]);
		v2_p=(v+DIM_N-s_minusone[i]);
		
		for(j=0;j<vec_num;j+=32)
		{
			tmp_sum1 = _mm256_load_si256((__m256i *)(sum1+j));
			tmp_v1   = _mm256_lddqu_si256((__m256i *)(v1_p+j));
			tmp_sum1 = _mm256_add_epi8(tmp_sum1, tmp_v1);
			_mm256_store_si256((__m256i*)(sum1+j),tmp_sum1);
			
			tmp_sum2 = _mm256_load_si256((__m256i *)(sum2+j));
			tmp_v2   = _mm256_lddqu_si256((__m256i *)(v2_p+j));
			tmp_sum2 = _mm256_add_epi8(tmp_sum2, tmp_v2);
			_mm256_store_si256((__m256i*)(sum2+j),tmp_sum2);
		}
	}
	
	for(i=0;i<vec_num;i++)
	{
		b[i]=sum1[i]-sum2[i];
	}
	
	return 0;
}


//b=as+e 
int poly_aff(const unsigned char *a, const unsigned char *s, unsigned char *e, unsigned char *b, unsigned int vec_num)
{
	int i;
	//b=ar
	poly_mul(a,s,b,vec_num);
	
	for(i=0;i<vec_num;i++)
	{
		b[i]+=e[i];
	}
	
	return 0;
}

// compress: cut the low 4bit
int poly_compress(const unsigned char *in,  unsigned char *out, const unsigned int vec_num)
{
	int i,loop;
	__m256i tmp0, tmp1;
	unsigned char *p_out;
	const unsigned char *p_in;
	
	if(vec_num>=64)
	{
		loop=vec_num-vec_num%64;
		for(i=0;i<loop;i+=64)
		{
			tmp0 = _mm256_lddqu_si256((__m256i *)(in+i));
			tmp1 = _mm256_lddqu_si256((__m256i *)(in+i+32));
			
			tmp0 = _mm256_srli_epi16(tmp0,4);
			tmp0 = _mm256_and_si256(tmp0,_mm256_set1_epi8(0x0F));
			tmp1 = _mm256_and_si256(tmp1,_mm256_set1_epi8(0xF0));
			tmp0 = _mm256_xor_si256(tmp1,tmp0);
			
			_mm256_storeu_si256((__m256i*)(out+(i>>1)),tmp0);
		}
	}
	
	loop=(vec_num%64)/2;
	p_in=in+vec_num-vec_num%64;
	p_out=out+(vec_num/64)*32;
	for(i=0;i<loop;i++)
	{
		p_out[i]=(p_in[i*2])>>4;
		p_out[i]=p_out[i]^((p_in[i*2+1])&0xf0);
	}
	
	return 0;
}
// decompress: set the low 4bits to be 0
int poly_decompress(const unsigned char *in,  unsigned char *out, const unsigned int vec_num)
{
	int i,loop;
	__m256i tmp0, tmp1;
	unsigned char *p_out;
	const unsigned char *p_in;
	
	if(vec_num>=64)
	{
		loop=vec_num-vec_num%64;
		for(i=0;i<loop;i+=64)
		{
			tmp0 = _mm256_lddqu_si256((__m256i *)(in+(i>>1)));
			
			tmp1 = _mm256_and_si256(tmp0,_mm256_set1_epi8(0xF0));
			tmp1 = _mm256_xor_si256(tmp1,_mm256_set1_epi8(0x08));
			
			tmp0 = _mm256_slli_epi16(tmp0,4);
			tmp0 = _mm256_and_si256(tmp0,_mm256_set1_epi8(0xF0));
			tmp0 = _mm256_xor_si256(tmp0,_mm256_set1_epi8(0x08));
			
			_mm256_storeu_si256((__m256i*)(out+i),tmp0);			
			_mm256_storeu_si256((__m256i*)(out+i+32),tmp1);
		}
	}
	
	loop=(vec_num%64)/2;
	p_in=in+(vec_num/64)*32;
	p_out=out+vec_num-vec_num%64;
	for(i=0;i<loop;i++)
	{
		p_out[i*2]=(p_in[i]<<4)^0x08;
		p_out[i*2+1]=(p_in[i]&0xf0)^0x08;
	}
	
	return 0;
}
