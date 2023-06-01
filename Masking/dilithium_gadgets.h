#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "random.h"
#include "convtable.h"
#include "convba_2014.h"
#include "utils.h"

#define uint128_t __uint128_t
#define int128_t __int128_t


void secMulAssignment(uint32_t* res, uint32_t* x, uint32_t q, int n);

void securediv(uint32_t* x, uint32_t* y, int log_alpha, int inv_alpha, int q, int n);


int ABC_rejection_sampling(uint32_t* x, int mode, int n);


void gen_y(uint32_t* y, int n);



void approximate_modulus_switching(uint64_t* x, uint32_t* y, int n);
void exact_modulus_switching(uint64_t* x, uint32_t* y, int n);



void decompose1(uint32_t* x, uint32_t* high, uint32_t* low, int n);
void decompose2(uint32_t* x, uint32_t* high, uint32_t* low, int n);