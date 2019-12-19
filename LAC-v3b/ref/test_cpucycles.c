#include <stdio.h>
#include <string.h>
#include <x86intrin.h>
#include "api.h"
#include "rand.h"
#include "ecc.h"
#include "bin-lwe.h"


#ifdef LINUX
#include <x86intrin.h>
#define cpucycles() _rdtsc()
#endif

#ifdef WIN
#include <intrin.h>
#define cpucycles() __rdtsc()
#endif

#define NTESTS 10000

unsigned long long sum,average;

static void print_uint64(unsigned long long num)
{
	if(num>=10)
		print_uint64(num/10);
	printf("%d",(int)(num%10));
}

int print_cpucycles(unsigned long long *t)
{
	int i;
	sum=0;
	for(i=0;i<NTESTS-1;i++)
	{
		t[i]=t[i+1]-t[i];
		sum+=t[i];
	}
	average=sum/(NTESTS-1);
	printf(": ");
	print_uint64(average);
    printf(" cpucycles\n");
	
	return 0;
}

//test pke
int test_pke_cpucycles()
{
	unsigned long long t[NTESTS];
	unsigned char pk[CRYPTO_PUBLICKEYBYTES];
	unsigned char sk[CRYPTO_SECRETKEYBYTES];
	unsigned char k1[CRYPTO_BYTES],k2[CRYPTO_BYTES],c[CRYPTO_CIPHERTEXTBYTES];
	int i;
	unsigned long long mlen=CRYPTO_BYTES,clen=CRYPTO_CIPHERTEXTBYTES;

    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_encrypt_keypair(pk,sk);
	}
	printf("pke kg     ");
	print_cpucycles(t);
	
	crypto_encrypt_keypair(pk,sk);
	random_bytes(k1,CRYPTO_BYTES);
	
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_encrypt(c,&clen,k1,mlen,pk);
	}
	
    printf("encryption ");
	print_cpucycles(t);
	
	crypto_encrypt_keypair(pk,sk);
	random_bytes(k1,CRYPTO_BYTES);
	crypto_encrypt(c,&clen,k1,mlen,pk);
	
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_encrypt_open(k2,&mlen,c,clen,sk);
	}
	
    printf("decryption ");
	print_cpucycles(t);
	printf("\n");
	
	return 0;
}

//test kem fo
int test_kem_fo_cpucycles()
{
	unsigned long long t[NTESTS];
	unsigned char pk[CRYPTO_PUBLICKEYBYTES];
	unsigned char sk[CRYPTO_SECRETKEYBYTES];
	unsigned char k1[CRYPTO_BYTES],k2[CRYPTO_BYTES],c[CRYPTO_CIPHERTEXTBYTES];
	int i;
	
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_kem_keypair(pk,sk);
	}
	printf("kem_fo kg  ");
	print_cpucycles(t);
	
	crypto_kem_keypair(pk,sk);
	random_bytes(k1,CRYPTO_BYTES);
	
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_kem_enc(c,k1,pk);
	}
	
    printf("encryption ");
	print_cpucycles(t);
	
	crypto_kem_keypair(pk,sk);
	crypto_kem_enc(c,k1,pk);
	
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_kem_dec(k2,c,sk);
	}
	
    printf("decryption ");
	print_cpucycles(t);
	printf("\n");
	
	return 0;
}

// test ke
int test_ke_cpucycles()
{
	unsigned long long t[NTESTS];
	unsigned char pk[CRYPTO_PUBLICKEYBYTES];
	unsigned char sk[CRYPTO_SECRETKEYBYTES];
	unsigned char k1[CRYPTO_BYTES],k2[CRYPTO_BYTES],c[CRYPTO_CIPHERTEXTBYTES];
	int i;
	
	//test alilce send
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_ke_alice_send(pk,sk);
	}
	
    printf("ke alice0  ");
	print_cpucycles(t);
	
	//test bob receive
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_ke_bob_receive(pk,c,k1);
	}
	
    printf("ke bob     ");
	print_cpucycles(t);
	
	//test alilce send
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_ke_alice_receive(pk,sk,c,k2);
	}
	
    printf("ke alice1  ");
	print_cpucycles(t);
	printf("\n");
	
	return 0;
}

// test ake
int test_ake_cpucycles()
{
	unsigned long long t[NTESTS];
	unsigned char pk_a[CRYPTO_PUBLICKEYBYTES],pk_b[CRYPTO_PUBLICKEYBYTES],pk[CRYPTO_PUBLICKEYBYTES];
	unsigned char sk[CRYPTO_SECRETKEYBYTES],sk_a[CRYPTO_SECRETKEYBYTES],sk_b[CRYPTO_SECRETKEYBYTES];
	unsigned char k_a[CRYPTO_BYTES],k_b[CRYPTO_BYTES],c_a[CRYPTO_CIPHERTEXTBYTES],c_b[2*CRYPTO_CIPHERTEXTBYTES],k1[CRYPTO_BYTES];
	int i;
	
	crypto_encrypt_keypair(pk_a,sk_a);
	crypto_encrypt_keypair(pk_b,sk_b);
	
	//test alilce send
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_ake_alice_send(pk,sk,pk_b,sk_a,c_a,k1);
	}
	
    printf("ake alice0 ");
	print_cpucycles(t);
	
	//test bob receive
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_ake_bob_receive(pk_b,sk_b,pk_a,pk,c_a,c_b,k_b);
	}
	
    printf("ake bob    ");
	print_cpucycles(t);
	
	//test alilce send
    for(i=0;i<NTESTS;i++)
	{
		t[i]=cpucycles();
		crypto_ake_alice_receive(pk_a,sk_a,pk_b,pk,sk,c_a,c_b,k1,k_a);
	}
	
    printf("ake alice1 ");
	print_cpucycles(t);
	printf("\n");

	return 0;
}

// test  hash
int test_hash_cpucycles()
{
	int i,loop=100000;
	unsigned char buf[DIM_N],seed[SEED_LEN],out[MESSAGE_LEN];
	uint64_t t0,t1;
	random_bytes(seed,SEED_LEN);
	pseudo_random_bytes(buf,DIM_N,seed);
	t0=__rdtsc();
	for(i=0;i<loop;i++)
	{
		hash(buf,DIM_N,out);
	}
	t1=__rdtsc();
	printf("test hash speed:\n");
//	printf("time: %f us\n",(t1-t0)/(2800.0*loop));
	printf("cpucycles: %f \n",(t1-t0)/((double)loop));
	printf("CPB: %f \n",(t1-t0)/((double)loop*DIM_N));
	printf("\n");
	
	return 0;
}

// test aes
int test_aes_cpucycles()
{
	int i,loop=100000;
	unsigned char buf[DIM_N],seed[SEED_LEN];
	uint64_t t0,t1;
	random_bytes(seed,SEED_LEN);
	t0=__rdtsc();
	for(i=0;i<loop;i++)
	{
		pseudo_random_bytes(buf,DIM_N,seed);
	}
	t1=__rdtsc();
	printf("test aes speed:\n");
//	printf("time: %f us\n",(t1-t0)/(2800.0*loop));
	printf("cpucycles: %f \n",(t1-t0)/((double)loop));
	printf("CPB: %f \n",(t1-t0)/((double)loop*DIM_N));
	printf("\n");
	
	return 0;
}

//test gen_psi
int test_gen_psi_cpucycles()
{
	int i,loop=10000;
	unsigned char seed[SEED_LEN];
	unsigned char e[DIM_N];
	uint64_t t0,t1;
	random_bytes(seed,SEED_LEN);
	
	t0=__rdtsc();
	for(i=0;i<loop;i++)
	{
		gen_e(e,seed);
	}
	t1=__rdtsc();
	printf("gen_psi_fix_ham cpucycles: %f \n",(t1-t0)/((double)loop));
	printf("\n");
	
	return 0;
}

//test gen_a
int test_gen_a_cpucycles()
{
	int i,loop=10000;
	unsigned char seed[SEED_LEN];
	unsigned char a[DIM_N];
	uint64_t t0,t1,sum;
	random_bytes(seed,SEED_LEN);
	sum=0;
	for(i=0;i<loop;i++)
	{
		t0=__rdtsc();
		gen_a(a,seed);
		t1=__rdtsc();
		sum+=(t1-t0);
	}
	
	printf("test gen_a speed:\n");
//	printf("time: %f us\n",(sum)/(2800.0*loop));
	printf("cpucycles: %f \n",(sum)/((double)loop));
	printf("\n");
	
	return 0;
}

//test polymul
int test_poly_mul_cpucycles()
{
	int i,loop=10000;
	unsigned char a[DIM_N],pk[DIM_N],seed[SEED_LEN];
	unsigned char sk[DIM_N];
	uint64_t t0,t1,sum;
	
	random_bytes(a,DIM_N);
	random_bytes(seed,SEED_LEN);
	gen_r(sk,seed);

	sum=0;
	for(i=0;i<loop;i++)
	{
		t0=__rdtsc();
		poly_mul(a,sk,pk,DIM_N);
		t1=__rdtsc();
		sum+=(t1-t0);
	}
	
	printf("test poly_nul speed:\n");
	printf("cpucycles: %f \n",(sum)/((double)loop));
	printf("\n");
	
	return 0;
}

int test_init()
{
	int i,loop=10000;
	unsigned char a[DIM_N],pk[DIM_N],seed[SEED_LEN];
	unsigned char sk[DIM_N];
	uint64_t t0,t1,sum;
	
	random_bytes(a,DIM_N);
	random_bytes(seed,SEED_LEN);
	gen_r(sk,seed);

	sum=0;
	for(i=0;i<loop;i++)
	{
		t0=__rdtsc();
		poly_mul(a,sk,pk,DIM_N);
		t1=__rdtsc();
		sum+=(t1-t0);
	}
	
	return 0;
}
