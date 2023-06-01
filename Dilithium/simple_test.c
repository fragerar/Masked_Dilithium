#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "params.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "symmetric.h"
#include "fips202.h"
#include "randombytes.h"
#include "sign.h"
#include "masked_sign.h"
#include "masking_interface.h"
#include "masked_polyvec_operations.h"
#include "packing.h"
#include "./test/cpucycles.h"

#define MLEN 59
#define NTESTS 10







int unmasked_test(void);
int masked_test(void);


int unmasked_test(void)
{
  unsigned int i, j;
  int ret;
  size_t mlen, smlen;
  uint8_t m[MLEN] = {0};
  uint8_t sm[MLEN + CRYPTO_BYTES];
  uint8_t m2[MLEN + CRYPTO_BYTES];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];

  for(i = 0; i < NTESTS; ++i) {
    randombytes(m, MLEN);

    crypto_sign_keypair(pk, sk);
    crypto_sign(sm, &smlen, m, MLEN, sk);
    ret = crypto_sign_open(m2, &mlen, sm, smlen, pk);

    if(ret) {
      fprintf(stderr, "Verification failed\n");
      return -1;
    }

    if(mlen != MLEN) {
      fprintf(stderr, "Message lengths don't match\n");
      return -1;
    }

    for(j = 0; j < mlen; ++j) {
      if(m[j] != m2[j]) {
        fprintf(stderr, "Messages don't match\n");
        return -1;
      }
    }

    randombytes((uint8_t *)&j, sizeof(j));
    do {
      randombytes(m2, 1);
    } while(!m2[0]);
    sm[j % CRYPTO_BYTES] += m2[0];
    ret = crypto_sign_open(m2, &mlen, sm, smlen, pk);
    if(!ret) {
      fprintf(stderr, "Trivial forgeries possible\n");
      return -1;
    }
  }
  printf("Unmasked Success !!\n");

  return 0;

}


int masked_test(void){
  unsigned int i, j;
  int ret;
  size_t mlen, smlen;
  uint8_t m[MLEN] = {0};
  uint8_t sm[MLEN + CRYPTO_BYTES];
  uint8_t m2[MLEN + CRYPTO_BYTES];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];

  const int ITER = 1000 ;
  int count = 0;

  for(i = 0; i < ITER; ++i) {
    printf("ITER: %i/%i\n", i+1, ITER);
    randombytes(m, MLEN);

    crypto_sign_keypair(pk, sk);



    uint8_t seedbuf[2*SEEDBYTES + 3*CRHBYTES];
    uint8_t *rho, *tr, *key;
    uint16_t nonce = 0;
    masked_polyvecl ms1;
    masked_polyveck ms2;
    polyvecl s1;
    polyveck s2, t0;

    rho = seedbuf;
    tr = rho + SEEDBYTES;
    key = tr + CRHBYTES;
    unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

    mask_polyvecl(&ms1, &s1);
    mask_polyveck(&ms2, &s2);


    count += masked_crypto_sign(sm, &smlen, m, MLEN, seedbuf, &ms1, &ms2, &t0);
    ret = crypto_sign_open(m2, &mlen, sm, smlen, pk);

    if(ret) {
      fprintf(stderr, "Verification failed\n");
      return -1;
    }

    if(mlen != MLEN) {
      fprintf(stderr, "Message lengths don't match\n");
      return -1;
    }

    for(j = 0; j < mlen; ++j) {
      if(m[j] != m2[j]) {
        fprintf(stderr, "Messages don't match\n");
        return -1;
      }
    }

    randombytes((uint8_t *)&j, sizeof(j));
    do {
      randombytes(m2, 1);
    } while(!m2[0]);
    sm[j % CRYPTO_BYTES] += m2[0];
    ret = crypto_sign_open(m2, &mlen, sm, smlen, pk);
    if(!ret) {
      fprintf(stderr, "Trivial forgeries possible\n");
      return -1;
    }
  }
  printf("Masked Success !!\n");

  printf("Avg repetitions: %f\n", ((double)count)/ITER);


  return 0;

}



void bench_unmasked(){
  const int ITER = 10000;
  unsigned int i, j;
  int ret;
  size_t mlen, smlen;
  uint8_t m[MLEN] = {0};
  uint8_t sm[MLEN + CRYPTO_BYTES];
  uint8_t m2[MLEN + CRYPTO_BYTES];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  crypto_sign_keypair(pk, sk);

  uint64_t global_start = cpucycles();

  for(i = 0; i < ITER; ++i) {
    randombytes(m, MLEN);
    crypto_sign(sm, &smlen, m, MLEN, sk);

  }
  uint64_t global_stop = cpucycles();


  printf("Total : %.0f \n", (double)(global_stop - global_start)/(ITER*1000));


}

void bench_signature_components(){

  unsigned int i;
  int ret;
  size_t mlen, smlen;
  uint8_t m[MLEN] = {0};
  uint8_t sm[MLEN + CRYPTO_BYTES];
  uint8_t m2[MLEN + CRYPTO_BYTES];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];


  const int ITER = 1000;

  uint8_t seedbuf[2*SEEDBYTES + 3*CRHBYTES];
  uint8_t *rho, *tr, *key;
  uint16_t nonce = 0;
  masked_polyvecl ms1;
  masked_polyveck ms2;
  polyvecl s1;
  polyveck s2, t0;

  crypto_sign_keypair(pk, sk);
  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + CRHBYTES;

  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
  randombytes(m, MLEN);

  uint64_t bench_vector[7] = {0,0,0,0,0,0,0};
  int count = 0;

  uint64_t global_start = cpucycles();

  for(int i=0; i < ITER; i++) {
    unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
    mask_polyvecl(&ms1, &s1);
    mask_polyveck(&ms2, &s2);

    randombytes(m, MLEN);
    count += bench_masked_crypto_sign(sm, &smlen, m, MLEN, seedbuf, &ms1, &ms2, &t0, bench_vector);


  }
  uint64_t global_stop = cpucycles();

  uint64_t total = 0;
  printf("Avg repetitions: %f\n", ((double)count)/ITER);
  for(int i=0; i < 7; ++i){
    printf("%.0f ", (double)bench_vector[i]/(ITER*1000)); 
  }
  printf("\n");

  printf("Total : %.0f \n", (double)(global_stop - global_start)/(ITER*1000));

}



int main(void){
  //unmasked_test();
  //masked_test();
  bench_signature_components();
  //bench_unmasked();

  return 0;
}
