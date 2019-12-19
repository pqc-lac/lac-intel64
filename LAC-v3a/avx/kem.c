#include <string.h>
#include "api.h"
#include "rand.h"
#include "bch.h"
#include "ecc.h"

//generate keypair
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
	//call the key generation algorithm of pke
	crypto_encrypt_keypair(pk, sk);
	return 0;
}
int crypto_kem_enc( unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{
	kem_enc_fo(pk,ss,ct);
	return 0;
}
int crypto_kem_dec( unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{
	const unsigned char *pk;
	pk=sk+SK_LEN;
	kem_dec_fo(pk,sk,ct,ss);
	return 0;
}
// fo encryption for cca security 
int kem_enc_fo(const unsigned char *pk, unsigned char *k, unsigned char *c)
{
	unsigned char buf[MESSAGE_LEN],seed[SEED_LEN],seed_buf[MESSAGE_LEN+SEED_LEN];
	unsigned long long clen;

	
	//check parameter
	if(pk==NULL || k==NULL || c==NULL)
	{
		return -1;
	}
	
	//generate random message m, stored in buf
	random_bytes(buf,MESSAGE_LEN);
	//compute seed=gen_seed(m|pk), add pk for multi key attack protection
	memcpy(seed_buf,buf,MESSAGE_LEN);
	memcpy(seed_buf+MESSAGE_LEN,pk,SEED_LEN);
	gen_seed(seed_buf,MESSAGE_LEN+SEED_LEN,seed);
	//encrypt m with seed
	pke_enc_seed(pk,buf,MESSAGE_LEN,c,&clen,seed);
	
	//compute k=hash(m|c)
	hash_to_k(buf,MESSAGE_LEN,k);
	
	return 0;
}

// fo encryption for cca security with seed
int kem_enc_fo_seed(const unsigned char *pk, unsigned char *k, unsigned char *c, unsigned char *seed)
{
	unsigned char buf[MESSAGE_LEN],local_seed[SEED_LEN],seed_buf[MESSAGE_LEN+SEED_LEN];
	unsigned long long clen;

	
	//check parameter
	if(pk==NULL || k==NULL || c==NULL)
	{
		return -1;
	}
	
	//generate random message m from seed, stored in buf
	pseudo_random_bytes(buf,MESSAGE_LEN,seed);
	//compute loacal_seed=gen_seed(m|pk), add pk for multi key attack protection
	memcpy(seed_buf,buf,MESSAGE_LEN);
	memcpy(seed_buf+MESSAGE_LEN,pk,SEED_LEN);
	gen_seed(seed_buf,MESSAGE_LEN+SEED_LEN,local_seed);
	//encrypt m with local_seed
	pke_enc_seed(pk,buf,MESSAGE_LEN,c,&clen,local_seed);
	
	//compute k=hash(m|c)
	hash_to_k(buf,MESSAGE_LEN,k);
	
	return 0;
}

// decrypt of fo mode
int kem_dec_fo(const unsigned char *pk, const unsigned char *sk, const unsigned char *c, unsigned char *k)
{
	unsigned char buf[MESSAGE_LEN+CIPHER_LEN],seed[SEED_LEN],seed_buf[MESSAGE_LEN+SEED_LEN];
	unsigned long long mlen,clen;
	unsigned char c_v[CIPHER_LEN];
	
	//check parameter
	if(pk==NULL || sk==NULL || k==NULL || c==NULL)
	{
		return -1;
	}
	
	//compute m from c
	pke_dec(sk,c,CIPHER_LEN, buf,&mlen);
	//compte k=hash(m|c)
	hash_to_k(buf,MESSAGE_LEN,k);
	//re-encryption with seed=gen_seed(m|pk), add pk for multi key attack protection
	memcpy(seed_buf,buf,MESSAGE_LEN);
	memcpy(seed_buf+MESSAGE_LEN,pk,SEED_LEN);
	gen_seed(seed_buf,MESSAGE_LEN+SEED_LEN,seed);
	pke_enc_seed(pk,buf,MESSAGE_LEN,c_v,&clen,seed);
	
	//verify
	if(memcmp(c,c_v,CIPHER_LEN)!=0)
	{
		//k=hash(hash(sk)|c)
		hash((unsigned char*)sk,SK_LEN,buf);
		hash(buf,MESSAGE_LEN+CIPHER_LEN,k);
	}
	
	return 0;
}

