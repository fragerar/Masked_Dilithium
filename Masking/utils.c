#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>

#include "random.h"
#include "utils.h"



int64_t cpucycles(void)
{ // Access system counter for benchmarking
  unsigned int hi, lo;

  asm volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
  return ((int64_t)lo) | (((int64_t)hi) << 32);
}


uint32_t genrand(int l)
{
  if (l==32) return rand32();
  return rand32() & ((1 << l)-1);
}


void share(uint32_t x,uint32_t a[],int n)
{
  a[0]=x;
  for(int i=1;i<n;i++)
    a[i]=0;
}

void share64(uint64_t x,uint64_t a[],int n)
{
  a[0]=x;
  for(int i=1;i<n;i++)
    a[i]=0;
}

void refreshBool(uint32_t a[],int l,int n)
{
  for(int i=0;i<n-1;i++)
  {
    uint32_t tmp=genrand(l);
    a[n-1]=a[n-1] ^ tmp;
    a[i]=a[i] ^ tmp;
  }
}


void refreshBool64(uint64_t a[],int n)
{
  for(int i=0;i<n-1;i++)
  {
    uint64_t tmp=rand64();
    a[n-1]=a[n-1] ^ tmp;
    a[i]=a[i] ^ tmp;
  }
}

void FullRefreshBool(uint32_t *a,uint32_t *b,int l,int n)
{
  for(int i=0;i<n;i++)
    b[i]=a[i];

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=genrand(l);
      b[i]^=tmp;
      b[j]^=tmp;
    }
  }
}

void refreshArith(uint32_t a[],int l,int n)
{
  uint32_t tmp;
  uint32_t ma=(1 << l)-1; 
  for(int i=0;i<n-1;i++)
  {
    tmp=rand32();
    a[n-1]=(a[n-1] + tmp) & ma;
    a[i]=(a[i] - tmp) & ma;
  }
}

void refreshArith64(uint64_t a[],int l,int n)
{
  uint64_t ma=((uint64_t)1 << l)-1; 
  for(int i=0;i<n-1;i++)
  {
    uint64_t tmp=rand64();
    a[n-1]=(a[n-1] + tmp) & ma;
    a[i]=(a[i] - tmp) & ma;
  }
}

void refreshArithModp(uint32_t a[],uint32_t p,int n)
{
  for(int i=0;i<n-1;i++)
  {
    uint32_t tmp=rand32() % p; //rand();
    a[n-1]=(a[n-1] + tmp) % p;
    a[i]=(a[i] + p - tmp) % p;
  }
}

void refreshArithModp64(uint64_t a[],uint64_t p,int n)
{
  for(int i=0;i<n-1;i++)
  {
    uint64_t tmp=rand64() % p; //rand();
    a[n-1]=(a[n-1] + tmp) % p;
    a[i]=(a[i] + p - tmp) % p;
  }
}

uint32_t xorop(uint32_t a[],int n)
{
  uint32_t r=0;
  for(int i=0;i<n;i++)
    r^=a[i];
  return r;
}

uint32_t addop(uint32_t a[],int l,int n)
{
  uint32_t r=0;
  for(int i=0;i<n;i++)
    r+=a[i];
  if (l==32)
    return r;
  else
    return r % (1 << l);
}

uint32_t addopmodp(uint32_t *a,int p,int n)
{
  uint32_t r=0;
  for(int i=0;i<n;i++)
    r+=a[i];
  return r % p;
}


void SecAnd(uint32_t *a,uint32_t *b,uint32_t *c,int k,int n)
{
  for(int i=0;i<n;i++)
    c[i]=a[i] & b[i];

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32(); //rand();
      uint32_t tmp2=(tmp ^ (a[i] & b[j])) ^ (a[j] & b[i]);
      c[i]^=tmp;
      c[j]^=tmp2;
    }
  }
  for(int i=0;i<n;i++) c[i]=c[i] % (1 << k);
}

void SecAdd(uint32_t *x,uint32_t *y,uint32_t *z,int k,int n)
{
  uint32_t u[n];
  for(int i=0;i<n;i++) u[i]=0;
  uint32_t w[n];
  SecAnd(x,y,w,k,n);
  uint32_t a[n];
  for(int i=0;i<n;i++) a[i]=x[i] ^ y[i];
  for(int j=0;j<k-1;j++)
  {
    uint32_t ua[n];
    SecAnd(u,a,ua,k,n);
    for(int i=0;i<n;i++) u[i]=(2*(ua[i] ^ w[i])) % (1 << k);
  }
  for(int i=0;i<n;i++) z[i]=x[i] ^ y[i] ^ u[i];
}



// Secure multiplication modulo 2^32
void SecMul(uint32_t *a,uint32_t *b,uint32_t *c,int n)
{
  for(int i=0;i<n;i++)
    c[i]=a[i]*b[i];

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32(); //rand();
      uint32_t tmp2=(tmp+a[i]*b[j])+a[j]*b[i];
      c[i]-=tmp;
      c[j]+=tmp2;
    }
  }
}

void SecMultModp(uint32_t *a,uint32_t *b,uint32_t *c,uint32_t p,int n)
{
  int i,j; 
  for(i=0;i<n;i++)
    c[i]=((uint64_t)a[i]*b[i]) % p;

  for(i=0;i<n;i++)
  {
    for(j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32() % p; //rand();
      uint32_t tmp2=((tmp+(uint64_t)a[i]*b[j])+(uint64_t)a[j]*b[i]) % p;
      c[i]=(c[i]+p-tmp) % p;
      c[j]=(c[j]+tmp2) % p;
    }
  }
}




void printShares(uint32_t *a,int n)
{
  for(int i=0;i<n;i++)
    printf("%u ",a[i]);
  printf("\n");
}

void printArith(uint32_t *a, int p, int n)
{
  for(int i=0;i<n;i++)
    printf("%u ",a[i]);
  if (p == 0) printf(" = %u \n", addop(a, 32, n));
  else printf(" = %u \n", addopmodp(a, p, n));
  
}

void printBool(uint32_t *a,int n)
{
  for(int i=0;i<n;i++)
    printf("%u ",a[i]);
  printf(" = %u \n", xorop(a, n));
}



void printArith64(uint64_t *a, uint64_t p, int n)
{
  uint64_t res = 0;
  for(int i=0;i<n;i++){
    if (p == 0) res = (res + a[i]);
    else        res = (res + a[i])%p;
    printf("%lu ",a[i]);
  }
  printf(" = %lu \n", res);
  
}

void printBool64(uint64_t *a,int n)
{
  uint64_t res = 0;
  for(int i=0;i<n;i++){
    res = (res ^ a[i]);
    printf("%lu ",a[i]);
  }
  printf(" = %lu \n", res);
}







