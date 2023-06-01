#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "random.h"
#include "convtable.h"
#include "convba_2014.h"
#include "utils.h"
#include "dilithium_gadgets.h"


#ifndef DILITHIUM_MODE
#define DILITHIUM_MODE 2
#endif


#define DIL_Q 8380417
#define DIL_Q_PRIME 134086672
#define ALPHA 4
//We assume alpha = 4 which means we can go up to 15 shares

#if DILITHIUM_MODE == 2
#define MU 18
#define K_APPROX 41
#define A_APPROX 262401
#define K_EXACT 45
#define A_EXACT 4198404
#define DELTA 44
#define D_BETA 78
#define D_GAMMA1 (1 << 17)


#elif DILITHIUM_MODE == 3
#define MU 20
#define K_APPROX 43
#define A_APPROX 1049601
#define K_EXACT 47
#define A_EXACT 16793615
#define DELTA 16
#define D_BETA 196
#define D_GAMMA1 (1 << 19)

#elif DILITHIUM_MODE == 5
#define MU 20
#define K_APPROX 43
#define A_APPROX 1049601
#define K_EXACT 47
#define A_EXACT 16793615
#define DELTA 16
#define D_BETA 120
#define D_GAMMA1 (1 << 19)
#endif


#ifdef COUNT
uint64_t count_rand = 0;
#endif



void print_shares_64(uint64_t* x, int p, int n){
  uint64_t sum = 0;
  for(int i=0; i < n; ++i) {
    printf("%lu ", x[i]);
    sum += x[i];
  }
  printf(" = %lu\n", sum%((uint64_t)1<<p));

}




void secMulAssignment(uint32_t* res, uint32_t* x, uint32_t q, int n){
  uint32_t c[n];

  for(int i=0;i<n;i++)
    c[i]= (res[i]*x[i])%q;

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32()%q;
      uint64_t tmp2=(tmp+res[i]*x[j])+res[j]*x[i];
      tmp2 %= q;
      c[i] = (c[i] + q - tmp)%q;
      c[j] = (c[j] +tmp2)%q;
    }
  }


  for(int i=0; i < n; ++i) res[i] = c[i];
}



void generic_1bit_shift(uint32_t* x, uint32_t* y, int q, int n){
  /* Shift of 1 bit from mod q to mod q/2 for any even q*/
  uint32_t b[n], a[n], z[n];

  for(int i=0; i < n; ++i) b[i] = x[i]&1;
  bool2ArithSPOGmodq(b, a, q, n);


  for(int i=0; i < n; ++i) z[i] = (x[i] + q - a[i])%q;

  for(int i=0; i < n-1; ++i){
    z[n-1] = (z[n-1]   + (z[i]&1))%q;
    z[i]   = (z[i] + q - (z[i]&1))%q;
  }
  for(int i=0; i < n; ++i) y[i] = z[i]>>1;
}

void generic_1bit_shift64(uint64_t* x, uint64_t* y, uint64_t q, int n){
  /* Shift of 1 bit from mod q to mod q/2 for any even q*/
  uint64_t b[n], a[n], z[n];

  for(int i=0; i < n; ++i) b[i] = x[i]&1;
  bool2ArithSPOGmodq64(b, a, q, n);


  for(int i=0; i < n; ++i) z[i] = (x[i] + q - a[i])%q;

  for(int i=0; i < n-1; ++i){
    z[n-1] = (z[n-1]   + (z[i]&1))%q;
    z[i]   = (z[i] + q - (z[i]&1))%q;
  }
  for(int i=0; i < n; ++i) y[i] = z[i]>>1;
}


void generic_shift(uint32_t* x, uint32_t* y, int k, int q, int n){
  /* Shift of k bits from mod 2^k * q to mod q for any q*/
  for(int i=0; i < (k>>1); ++i){
    generic_1bit_shift(x, y, (1<<(k-2*i))*q, n);
    generic_1bit_shift(y, x, (1<<(k-(2*i+1)))*q, n);
  }
  if (k&1){
    generic_1bit_shift(x, y, 2*q, n);
  } else {
    for(int i=0; i < n; ++i) y[i] = x[i];
  }

}

void generic_shift64(uint64_t* x, uint64_t* y, int k, uint64_t q, int n){
  /* Shift of k bits from mod 2^k * q to mod q for any q*/
  for(int i=0; i < (k>>1); ++i){
    generic_1bit_shift64(x, y, (1<<(k-2*i))*q, n);
    generic_1bit_shift64(y, x, (1<<(k-(2*i+1)))*q, n);
  }
  if (k&1){
    generic_1bit_shift64(x, y, 2*q, n);
  } else {
    for(int i=0; i < n; ++i) y[i] = x[i];
  }

}



int ABC_rejection_sampling(uint32_t* x, int mode, int n){
  /* Assumes input masked mod+ q, returns 1 if  value is valid*/


  const int BITSIZE = 24;
  const int beta = D_BETA;
  uint32_t a=0, res = 0;
  if (mode == 0)  a = D_GAMMA1 - beta;
  else if (mode == 1) a = ((DIL_Q-1)/(2*DELTA)) - beta;

  uint32_t upper[n], b[n], y[n];
  for(int i=1; i < n; ++i) upper[i] = 0;
  upper[0] = (~((a<<1)-1)+1)%(1<<BITSIZE);


  x[0] = (x[0] + a - 1)%DIL_Q;


  ConvertABModp(x, b, DIL_Q, BITSIZE, n);

  SecAdd(b, upper, y, BITSIZE, n);

  for(int i=0; i < n; ++i) res ^= (y[i] >> (BITSIZE-1));

  return res;


}


void rejection_sampling(uint32_t* x, uint32_t* res, int mode, int n){
  const int beta = 196;
  uint32_t a=0;
  if (mode == 0)      a = (1<<17) - beta;
  else if (mode == 1) a = ((DIL_Q-1)/32) - beta;
  const int p = 10687;
  const int l = 10; //p*2^l = 10943488 > q
  int rho = l + ALPHA; // 14 (p*2^rho is 28 bits)
  int modulus = p*(1<<rho);
  int u = ((int64_t)modulus*a/DIL_Q);


  uint32_t y[n], ny[n], positive_side[n], negative_side[n];
  y[0] = ((((int64_t)x[0]*modulus)/DIL_Q)+n-1)%modulus;
  for(int i=1; i < n; ++i) y[i] = (((int64_t)x[i]*modulus)/DIL_Q)%modulus;

  ny[0] = ((((int64_t)(DIL_Q-x[0])*modulus)/DIL_Q)+n-1)%modulus;
  for(int i=1; i < n; ++i) ny[i] = (((int64_t)(DIL_Q-x[i])*modulus)/DIL_Q)%modulus;


  y[0] = (y[0] + modulus - u)%modulus;


  generic_shift(y, positive_side, rho, p, n);



  ny[0] = (ny[0] + modulus - u)%modulus;
  generic_shift(ny, negative_side, rho, p, n);


  SecMultModp(positive_side, negative_side, res, p, n);

}



void gen_y(uint32_t* y, int n){
  int k=K_EXACT;
  uint64_t x[n];
  uint64_t arith_x[n];

  for(int i=0; i < n; ++i) x[i] = rand32()%(1<<MU);
  impconvBA64(arith_x, x, n);
  for(int i=0; i < n; ++i) arith_x[i] %= (1<<k);


  exact_modulus_switching(arith_x, y, n);
  y[0] = (y[0] + DIL_Q - (1<<(MU-1)))%DIL_Q;

}


void approximate_modulus_switching(uint64_t* x, uint32_t* y, int n){

  uint32_t k = K_APPROX;
  uint32_t a = A_APPROX;

  uint128_t temp;
  

  y[0] = (((((uint128_t)x[0]*a*DIL_Q))>>k) + n - 1)%DIL_Q;
  for(int i=1; i < n; ++i){
    temp  = (uint128_t)x[i]*a*DIL_Q;
    temp >>= k;
    y[i] = temp % DIL_Q;
  }
}


void exact_modulus_switching(uint64_t* x, uint32_t* y, int n){
  //assume alpha = 4
  uint32_t k = K_EXACT;
  uint32_t a = A_EXACT;
  uint128_t temp;
  uint32_t z[n];
  

  for(int i=0; i < n; ++i){
    temp  = (uint128_t)x[i]*a*DIL_Q;
    temp >>= (k-ALPHA);
    z[i] = (uint32_t)(temp % DIL_Q_PRIME);
  }
  z[0] = (z[0] + n - 1)%DIL_Q_PRIME;
  generic_shift(z, y, ALPHA, DIL_Q, n);
}



void AB_convert_shift3(uint32_t* x, uint32_t* y, int k, int n){
  //Arithmetic to Boolean conversion of k bits.
  uint32_t temp[n];
  for(int i=0; i < n; ++i) y[i] = 0;

  for (int i=0; i < k; ++i){
    for(int j=0; j < n; ++j) y[j] |= (x[j]&1)<<i;
    shift3(x, temp, k, n);
    for(int j=0; j < n; ++j) x[j] = temp[j];
  }

}



void decompose1(uint32_t* x, uint32_t* high, uint32_t* low, int n){
  const int RHO = 25; //assume n <= 8.

  

  #if ((DILITHIUM_MODE == 3) || (DILITHIUM_MODE == 5))
  uint32_t modulus = (DELTA << RHO);
  uint32_t y[n], z[n];


  for(int i=0; i < n; ++i) y[i] = ((((uint64_t)x[i]*DELTA)<<RHO)/DIL_Q)%modulus;
  y[0] = (y[0] + n - 1 + (1<<(RHO-1)))%modulus;

  generic_shift(y, z, RHO, DELTA, n);

  refreshArithModp(z, DELTA, n);

  #else

  uint64_t modulus = (DELTA << RHO);
  uint64_t y[n], z[n];

  for(int i=0; i < n; ++i) y[i] = ((((uint64_t)x[i]*DELTA)<<RHO)/DIL_Q)%modulus;
  y[0] = (y[0] + n - 1 + (1<<(RHO-1)))%modulus;
  generic_shift64(y, z, RHO, DELTA, n);
  refreshArithModp64(z, DELTA, n);

  #endif


  *high = 0;
  for(int i=0; i < n; ++i) *high = (*high + z[i])%DELTA;
  
  for(int i=0; i < n; ++i) low[i] = x[i];
  low[0] = (low[0] + DIL_Q - *high*(DIL_Q-1)/DELTA)%DIL_Q;
  

}


void decompose2(uint32_t* x, uint32_t* high, uint32_t* low, int n){

  const int BITSIZE = 24;

  uint32_t s[n], b[n];
  #if (DILITHIUM_MODE == 2)
  uint32_t a[n], t[n], masked_high[n];
  #endif




  for(int i=0; i < n; ++i) s[i] = ((DIL_Q-DELTA)*(uint64_t)x[i])%DIL_Q;
  s[0] = (s[0] + (DIL_Q-1)/2)%DIL_Q;


  
  ConvertABModp(s, b, DIL_Q, BITSIZE, n);


  #if (DILITHIUM_MODE == 2)
  for(int j=0; j < n; ++j) t[j] = b[j]&1;
  bool2ArithSPOGmodq(t, masked_high, DELTA, n);
  
  for(int i=1; i < BITSIZE; ++i){
    for(int j=0; j < n; ++j) t[j] = (b[j] >> i)&1;
    bool2ArithSPOGmodq(t, a, DELTA, n);
    for(int j=0; j < n; ++j) masked_high[j] = (masked_high[j] + a[j] * (1<<i)%DELTA)%DELTA;
  }


  *high = 0;
  for(int i=0; i < n; ++i) *high = (*high + masked_high[i])%DELTA;
  #else
  *high = 0;
  for(int i=0; i < n; ++i) *high = (*high ^ (b[i]&0xF));
  #endif


  for(int i=0; i < n; ++i) low[i] = x[i];
  low[0] = (low[0] + DIL_Q - *high*(DIL_Q-1)/DELTA)%DIL_Q;
 
  
  
}

#define ITER 1000

void randomness_usage(){
  uint32_t max_n = 8;
  int64_t val;
  uint32_t x32[max_n];
  uint64_t x64[max_n];
  uint32_t y[max_n];

  uint32_t high, low[max_n];


  printf("Randomness CGV14: ");
  for(int j=2; j < max_n; ++j){
    val = rand();
    share(val, x32, j);
    refreshArith(x32, 32, j);
    count_rand = 0;
    ConvertAB(x32,y,32, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");
  


  printf("Randomness AB_convert_shift3: ");
  for(int j=2; j < max_n; ++j){
    val = rand();
    share((uint32_t)val, x32, j);
    refreshArith(x32, 32, j);
    count_rand = 0;
    AB_convert_shift3(x32, y, 32, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");





  printf("Randomness decompose1: ");
  for(int j=2; j < max_n; ++j){
    val = rand()%DIL_Q;
    share((uint32_t)val, x32, j);
    refreshArithModp(x32, DIL_Q, j);
    count_rand = 0;
    decompose1(x32, &high, low, j);
    printf("$%lu$ &", count_rand);
  }
  printf("\n");

  printf("Randomness decompose2: ");
  for(int j=2; j < max_n; ++j){
    val = rand()%DIL_Q;
    share((uint32_t)val, x32, j);
    refreshArithModp(x32, DIL_Q, j);
    count_rand = 0;
    decompose2(x32, &high, low, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");


  val = rand32();

  printf("Randomness SPOG19: ");
  for(int j=2; j < max_n; ++j){
    share(val, x32, j);
    refreshBool(x32, MU, j);
    count_rand = 0;
    bool2ArithSPOGmodqMulti(x32, y, MU, DIL_Q, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");


  
  printf("Randomness approx_mod_switch: ");
  for(int j=2; j < max_n; ++j){
    share64((uint64_t)val, x64, j);
    refreshArith64(x64, K_APPROX, j);
    count_rand = 0;
    approximate_modulus_switching(x64, y, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");


  printf("Randomness exact_mod_switch: ");
  for(int j=2; j < max_n; ++j){
    share64((uint64_t)val, x64, j);
    refreshArith64(x64, K_EXACT, j);
    count_rand = 0;
    exact_modulus_switching(x64, y, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");


  printf("Randomness impconv_BA: ");
  for(int j=2; j < max_n; ++j){
    share(val, x32, j);
    count_rand = 0;
    impconvBA(x32, y, j);
    printf(" $%lu$ &", count_rand);
  }
  printf("\n");





}


void bench_decompose(){

  uint64_t start, stop;
  uint32_t max_n = 8;
  int64_t val;
  uint32_t x32[max_n];


 uint32_t high, low[max_n];

  printf("Avg speed decomposeComp: ");
  for(int j=2; j < max_n; ++j){
    val = rand()%DIL_Q;
    share((uint32_t)val, x32, j);
    refreshArithModp(x32, DIL_Q, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) decompose1(x32, &high, low, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");

  printf("Avg speed decomposeMod: ");
  for(int j=2; j < max_n; ++j){
    val = rand()%DIL_Q;
    share((uint32_t)val, x32, j);
    refreshArithModp(x32, DIL_Q, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) decompose2(x32, &high, low, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");

}


void bench_conv_AB(){


  int BITSIZE = 32;
  uint64_t start, stop;
  uint32_t max_n = 11;
  int32_t val;
  uint32_t x32[max_n];
  uint32_t y[max_n];

  val = rand32();
  
  printf("Avg CGV14: ");
  for(int j=2; j < max_n; ++j){
    share(val, x32, j);
    refreshArith(x32, BITSIZE, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) ConvertAB(x32,y,BITSIZE, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");

  printf("Avg AB_convert_shiftmod: ");
  for(int j=2; j < max_n; ++j){
    share(val, x32, j);
    refreshArith(x32, BITSIZE, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) AB_convert_shift3(x32, y, BITSIZE, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


}



void bench_mu_bits_to_modq(){
  uint64_t start, stop;
  int max_n = 8;
  uint32_t x[max_n], y[max_n];
  uint64_t x64[max_n];
  uint32_t val;
  

  val = rand()%(1<<MU);


  printf("Avg speed SPOG19: ");
  for(int j=2; j < max_n; ++j){
    share(val, x, j);
    refreshBool(x, MU, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++)bool2ArithSPOGmodqMulti(x, y, MU, DIL_Q, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  
  printf("Avg speed approx_mod_switch: ");
  for(int j=2; j < max_n; ++j){
    share64((uint64_t)val, x64, j);
    refreshArith64(x64, K_APPROX, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) approximate_modulus_switching(x64, y, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  printf("Avg speed exact_mod_switch: ");
  for(int j=2; j < max_n; ++j){
    share64((uint64_t)val, x64, j);
    refreshArith64(x64, K_EXACT, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++) exact_modulus_switching(x64, y, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");


  printf("Avg speed impconv_BA: ");
  for(int j=2; j < max_n; ++j){
    share(val, x, j);
    start = cpucycles();
    for(int i=0; i < ITER; i++)  impconvBA(x, y, j);
    stop = cpucycles();
    //printf(" %f (%i) |", (double)(stop-start)/(ITER), j);
    printf(" $%.0f$ &", (double)(stop-start)/(ITER));
  }
  printf("\n");

}


#ifdef TESTS_SAMPLE_Y

int test_approx_mod_switch(){
  int n = 4;
  int k = K_APPROX;
  uint64_t x[n];
  uint32_t y[n];
  unsigned val, res;

  printf("Test approx_mod_switch...");
  for(int i=0; i < ITER; ++i){
    val = rand32()%(1<<MU);
    share64((uint64_t)val, x, n);    
    refreshArith64(x, k, n);
    approximate_modulus_switching(x, y, n);
    res = 0;
    for(int j=0; j < n; ++j) res = (res + y[j])%DIL_Q;
    if ((res < val) || res > val+(1<<ALPHA)) {
      printf("Fail! Iteration %i: (val, res) = (%i, %i)\n",i, val, res);
      return 0;
    }
  }
  printf("Success!\n");
  return 1;
}

int test_exact_mod_switch(){
  int n = 4;
  int k = K_EXACT;
  uint64_t x[n];
  uint32_t y[n];
  unsigned val, res;

  printf("Test exact_mod_switch...");
  for(int i=0; i < ITER; ++i){
    val = rand32()%(1<<MU);
    share64((uint64_t)val, x, n);    
    refreshArith64(x, k, n);
    exact_modulus_switching(x, y, n);
    res = 0;
    for(int j=0; j < n; ++j) res = (res + y[j])%DIL_Q;
    if (res != val) {
      printf("Fail! Iteration %i: (val, res) = (%i, %i)\n",i, val, res);
      return 0;
    }
  }
  printf("Success!\n");
  return 1;
}



int test_AB_convert_shift3(){
  int n = 4;
  int k = 30;
  uint32_t x[n];
  uint32_t y[n];
  unsigned val, res;

  printf("Test AB_convert_shift3...");
  for(int i=0; i < ITER; ++i){
    val = rand32()&((1LLU<<k)-1);
    share(val, x, n);    
    refreshArith(x, k, n);
    AB_convert_shift3(x, y, k, n);
    res = xorop(y, n);
    if (res != val) {
      printf("Fail! Iteration %i: (val, res) = (%i, %i)\n",i, val, res);
      return 0;
    }
  }
  printf("Success!\n");
  return 1;

}



int test_RS_z(){
  int n=6;
  const int beta = D_BETA;
  const int b = D_GAMMA1 + beta;
  const int bound = D_GAMMA1 - beta;
  //const int p = 10687;
  uint32_t x[n]; 
  int t;

  printf("Test RejectZ...");

  for(int val=0; val < b+1; val++){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    t = ABC_rejection_sampling(x, 0, n);

    if ((t == 0) && (val < bound)){
      printf("Fail. Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val >= bound)){
      printf("Fail. Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  } 



  for(int val = DIL_Q-1; val > DIL_Q - b - 1; --val){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    t = ABC_rejection_sampling(x, 0, n);
    if ((t == 0) && (val > DIL_Q - bound)){
      printf("Fail. Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val <= DIL_Q - bound)){
      printf("Fail. Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  }

  printf("Success\n");

  return 1;

}


int test_RS_r(){
  int n=6;
  const int beta = D_BETA;
  const int b = (DIL_Q-1)/(2*DELTA) + beta;
  const int bound = (DIL_Q-1)/(2*DELTA) - beta;
  uint32_t x[n];
  int t;

  printf("Test RejectR...");
  for(int val=0; val < b+1; val++){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);

    t = ABC_rejection_sampling(x, 1, n);
    if ((t == 0) && (val < bound)){
      printf("Fail. Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val >= bound)){
      printf("Fail. Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  } 

 

  for(int val = DIL_Q-1; val > DIL_Q - b - 1; --val){
    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    t = ABC_rejection_sampling(x, 1, n);
    if ((t == 0) && (val > DIL_Q - bound)){
      printf("Fail. Valid value rejected. Value: %i\n", val);
      return 0;
    }
    if ((t != 0) && (val <= DIL_Q - bound)){
      printf("Fail. Invalid value accepted. Value: %i\n", val);
      return 0;
    }
  }

  printf("Success\n");

  return 1;

}


int test_generic_shift(){
  int n = 6;
  uint32_t x[n], y[n], val;
  int k; 
  int q;
  uint32_t res;


  for(int j=0; j < ITER; ++j){
    res = 0;
    q = (rand()%49999) + 1;
    k = (rand()%9) + 1;

    val = rand()%((1<<k)*q);
    share(val, x, n);
    refreshArithModp(x, (1<<k)*q, n);
    generic_shift(x, y, k, q, n);

    for(int i=0; i < n; ++i) res = (res + y[i])%q;
    if (res != (val>>k)) {
      printf("Fail! Iteration %i: (val, res) = (%i, %i)\n",j, val>>k, res);
      printArith(x, (1<<k)*q, n);
      printArith(y, q, n);
      return 0;
    }
  }
  printf("Success!\n");



  return 1;
}




static int32_t dilithium_decompose(int32_t *a0, int32_t a) {
  int GAMMA2_TESTS = ((DIL_Q-1)/(2*DELTA));


  int32_t a1;

  a1  = (a + 127) >> 7;
#if DELTA == 16
  a1  = (a1*1025 + (1 << 21)) >> 22;
  a1 &= 15;
#elif DELTA == 44
  a1  = (a1*11275 + (1 << 23)) >> 24;
  a1 ^= ((43 - a1) >> 31) & a1;
#endif 

  *a0  = a - a1*2*GAMMA2_TESTS;
  *a0 -= (((DIL_Q-1)/2 - *a0) >> 31) & DIL_Q;
  return a1;
}


int test_decompose(){

  int n=5;
  int32_t val;
  int32_t r1, r0;
  uint32_t x[n], res_low[n], test;
  uint32_t res_high;



  printf("Test Decompose...");
  for(int i=0; i < DIL_Q; ++i){
    val = i;


    share(val, x, n);
    refreshArithModp(x, DIL_Q, n);
    
    decompose2(x, &res_high, res_low, n);

    r1 = dilithium_decompose(&r0, val);
    test = addopmodp(res_low, DIL_Q, n);


    if ((res_high != r1) || (test != (r0+DIL_Q)%DIL_Q)){
      printf("Fail!\n");
      return 0;
    }
  }
  printf("Success.\n");
  
  return 1;

}




void run_tests(){
  test_approx_mod_switch();
  test_exact_mod_switch();
  test_AB_convert_shift3();
  test_decompose();
  test_RS_z();
  test_RS_r();
} 


int main(){
  
  run_tests();
  
  bench_decompose();
  bench_mu_bits_to_modq();
  bench_conv_AB();
  randomness_usage();

}
#endif


