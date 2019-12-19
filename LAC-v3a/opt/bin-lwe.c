#include <string.h>
#include "bin-lwe.h"
#include "rand.h"
#include "lac_param.h"
#include "stdio.h"

//generate the public parameter a from seed
int gen_a(unsigned char *a,  const unsigned char *seed)
{
	int i,j;
	unsigned char buf[SEED_LEN];
	//check the pointers
	if(a==NULL || seed==NULL)
	{
		return 1;
	}
	
	pseudo_random_bytes(a,DIM_N,seed);
	
	hash(seed,SEED_LEN,buf);
	j=0;
	for(i=0;i<DIM_N;i++)
	{
		while(a[i]>=Q)
		{
			
			memcpy(a+i,buf+(j++),1);//replace a[i] with buf[j]
			if(j>=MESSAGE_LEN)
			{
				hash(buf,MESSAGE_LEN,buf);//use hash chain to refresh buf
				j=0;
			}
			
		}
	}
	
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

static int gen_index(uint16_t *index, unsigned char *seed)
{
	int i,rbuf_len;
	rbuf_len=SAMPLE_LEN+SAMPLE_LEN*(DIM_N>>8)/16;
	unsigned char  rbuf[DIM_N*2],*p;
	
	pseudo_random_bytes(rbuf,rbuf_len,seed);
	p=rbuf+SAMPLE_LEN;
	
	#if defined(LAC128) || defined(LAC_LIGHT)
	//prepare random index
	for(i=0;i<SAMPLE_LEN;i+=8)
	{
		index[i]=(rbuf[i]^(((*p)<<1)&0x100));
		index[i+1]=(rbuf[i+1]^(((*p)<<2)&0x100));
		index[i+2]=(rbuf[i+2]^(((*p)<<3)&0x100));
		index[i+3]=(rbuf[i+3]^(((*p)<<4)&0x100));
		index[i+4]=(rbuf[i+4]^(((*p)<<5)&0x100));
		index[i+5]=(rbuf[i+5]^(((*p)<<6)&0x100));
		index[i+6]=(rbuf[i+6]^(((*p)<<7)&0x100));
		index[i+7]=(rbuf[i+7]^(((*p)<<8)&0x100));
		p++;
	}
	#else
	for(i=0;i<SAMPLE_LEN;i+=4)
	{
		index[i]=(rbuf[i]^(((*p)<<2)&0x300));
		index[i+1]=(rbuf[i+1]^(((*p)<<4)&0x300));
		index[i+2]=(rbuf[i+2]^(((*p)<<6)&0x300));
		index[i+3]=(rbuf[i+3]^(((*p)<<2)&0x300));
		p++;
	}
	#endif
	
	return 0;
}

//for r and s
int gen_r(unsigned char *r,  unsigned char *seed)
{
	
	if(r==NULL)
	{
		return 1;
	}
	
	int i,p;
	uint16_t tmp;
	uint16_t  r_buf[DIM_N],index[SAMPLE_LEN],tmp_index,index_mk;
	unsigned int mask_p,loop=SAMPLE_LEN;

	//init r to be 1,2,3,4,5
	for(i=0;i<DIM_N;i++)
	{
		r_buf[i]=i;
	}
	
	
	p=0;
	while(p<NUM_ONE*2)
	{
		gen_index(index,seed);
	    //shuffle
		for(i=0;i<loop;i++)
		{
			//check index
			index_mk=index[i];
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
int mul_core(const unsigned char *a, const unsigned char *s, int64_t *sum1, int64_t *sum2, unsigned int vec_num)
{
	int i,j,loop;
	int16_t v[DIM_N+DIM_N];
	int64_t *v1_p,*v2_p;
	uint16_t *s_one,*s_minusone;
	
	//construct matrix of a
	for(i=0;i<DIM_N;i++)
	{
		v[i]=Q-a[i];
		v[i+DIM_N]=a[i];
	}	
	
	s_one=(uint16_t*)s;
	s_minusone=(uint16_t*)s+NUM_ONE;
	memset(sum1,0,vec_num*sizeof(int64_t));
	memset(sum2,0,vec_num*sizeof(int64_t));
	
	loop=vec_num/4;
	for(i=0;i<NUM_ONE;i++)
	{
		v1_p=(int64_t*)(v+DIM_N-s_one[i]);
		for(j=0;j<loop;j++)
		{
			sum1[j]+=v1_p[j];
		}
		
		v2_p=(int64_t*)(v+DIM_N-s_minusone[i]);
		for(j=0;j<loop;j++)
		{
			sum2[j]+=v2_p[j];
		}
	}
	
	return 0;
}
 

int poly_mul(const unsigned char *a, const unsigned char *s, unsigned char *b, unsigned int vec_num)
{
	int i,loop;
	int64_t sum1[DIM_N],sum2[DIM_N];
	
	mul_core(a,s,sum1,sum2,vec_num);
	
	loop=vec_num/4;
	for(i=0;i<loop;i++)
	{
		b[i*4  ]=((( sum1[i]     &0xffff)-( sum2[i]     &0xffff)+BIG_Q)%Q);
		b[i*4+1]=((((sum1[i]>>16)&0xffff)-((sum2[i]>>16)&0xffff)+BIG_Q)%Q);
		b[i*4+2]=((((sum1[i]>>32)&0xffff)-((sum2[i]>>32)&0xffff)+BIG_Q)%Q);
		b[i*4+3]=((((sum1[i]>>48)&0xffff)-((sum2[i]>>48)&0xffff)+BIG_Q)%Q);
	}
	
	return 0;
}
//b=as+e 
int poly_aff(const unsigned char *a, const unsigned char *s, unsigned char *e, unsigned char *b, unsigned int vec_num)
{
	int i,loop;
	int64_t sum1[DIM_N],sum2[DIM_N];
	
	mul_core(a,s,sum1,sum2,vec_num);
	
	loop=vec_num/4;

	for(i=0;i<loop;i++)
	{
		b[i*4  ]=((( sum1[i]     &0xffff)-( sum2[i]     &0xffff)+e[4*i]  +BIG_Q)%Q);
		b[i*4+1]=((((sum1[i]>>16)&0xffff)-((sum2[i]>>16)&0xffff)+e[4*i+1]+BIG_Q)%Q);
		b[i*4+2]=((((sum1[i]>>32)&0xffff)-((sum2[i]>>32)&0xffff)+e[4*i+2]+BIG_Q)%Q);
		b[i*4+3]=((((sum1[i]>>48)&0xffff)-((sum2[i]>>48)&0xffff)+e[4*i+3]+BIG_Q)%Q);
	}
	
	return 0;
}
 

// compress: cut the low 4bit
int poly_compress(const unsigned char *in,  unsigned char *out, const unsigned int vec_num)
{
	int i,loop;
	uint64_t tmp0, tmp1,*p64_in,*p64_out;
	unsigned char *p_out;
	const unsigned char *p_in;
	
	if(vec_num>=16)
	{
		loop=vec_num/8;
		p64_in=(uint64_t*)in;
		p64_out=(uint64_t*)out;
		for(i=0;i<loop;i+=2)
		{
			tmp0 = p64_in[i];
			tmp1 = p64_in[i+1];
			
			tmp0 >>= 4;
			tmp0 &= 0x0F0F0F0F0F0F0F0F;
			tmp1 &= 0xF0F0F0F0F0F0F0F0;
			tmp0 ^= tmp1;
			
			p64_out[i>>1]=tmp0;
		}
	}
	
	loop=(vec_num%16)/2;
	p_in=in+vec_num-vec_num%16;
	p_out=out+(vec_num/16)*8;
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
	uint64_t tmp0, tmp1,*p64_in,*p64_out;
	unsigned char *p_out;
	const unsigned char *p_in;
	
	if(vec_num>=16)
	{
		loop=vec_num/8;
		p64_in=(uint64_t*)in;
		p64_out=(uint64_t*)out;
		for(i=0;i<loop;i+=2)
		{
			tmp0 = p64_in[i>>1];
			
			tmp1 = (tmp0&0xF0F0F0F0F0F0F0F0)^0x0808080808080808;
						
			tmp0 <<= 4;
			tmp0 &= 0xF0F0F0F0F0F0F0F0;
			tmp0 ^= 0x0808080808080808;
			
			p64_out[i]=tmp0;
			p64_out[i+1]=tmp1;
		}
	}
	
	loop=(vec_num%16)/2;
	p_in=in+(vec_num/16)*8;
	p_out=out+vec_num-vec_num%16;
	for(i=0;i<loop;i++)
	{
		p_out[i*2]=(p_in[i]<<4)^0x08;
		p_out[i*2+1]=(p_in[i]&0xf0)^0x08;
	}
	
	return 0;
}



