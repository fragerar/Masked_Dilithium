// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.

// The file contains:
// 1) Arithmetic to Boolean conversion modulo 2^k from [CGV14]
// 2) Boolean to arithmetic conversion modulo 2^k from [CGV14]
// 3) Arithmetic to Boolean conversion modulo p, using the extended approach from [BBE+18]
// 4) Boolean to arithmetic conversion modulo p, using the extended approach from [BBE+18]
// 5) Boolean to arithmetic conversion modulo p from [SPOG19]

// For these algorithms, we use n=t+1 shares, while 1) and 2) where designed for n=2t+1 shares
// One should see in which case a refresh is necessary.

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include "utils.h"
#include "random.h"





uint32_t AddGoubin(uint32_t x,uint32_t y,int k)
{
  uint32_t u=0;
  uint32_t w= x & y;
  uint32_t a= x ^ y;
  for(int i=0;i<k-1;i++)
    u=2*( (u & a) ^ w);
  return a ^ u;
}

// We must have p<2^31
uint32_t AddGoubinModp(uint32_t x,uint32_t y,int k,uint32_t p)
{
  uint32_t m=-1;
  uint32_t s=AddGoubin(x,y,k);
  uint32_t sp=AddGoubin(s,-p,k);
  uint32_t b=-(sp >> (k-1));
  return (s & b) ^ (sp & (m ^ b)); 
}



void ExpandArith(uint32_t *x,uint32_t *xp,int n2,int n)
{
  for(int i=0;i<n/2;i++)
  {
    uint32_t r=rand32();
    xp[2*i]=x[i]-r;
    xp[2*i+1]=r;
  }
  if ((n & 1)==1) 
  {
    if (n2==n/2)
      xp[n-1]=0;
    else
      xp[n-1]=x[n2-1];
  }
}


// We assume that p<2^(k-1)
void SecAddModp(uint32_t *x,uint32_t *y,uint32_t *z,uint32_t p,int k,int n)
{
  uint32_t s[n];
  SecAdd(x,y,s,k,n);  // s=x+y

  uint32_t mp[n];
  share((1 << k)-p,mp,n);

  uint32_t sp[n];
  SecAdd(s,mp,sp,k,n);  // sp=x+y-p

  uint32_t b[n];
  for(int i=0;i<n;i++)
    b[i]=-(sp[i] >> (k-1));  // b=1...1 if sp<0. b=0 if sp>=0
  for(int i=0;i<n;i++) b[i]=b[i] % (1 << k);

  uint32_t c[n];
  FullRefreshBool(b,c,k,n);

  uint32_t z2[n];
  SecAnd(s,c,z2,k,n);    // if sp<0, then b=1....1 and we select s

  uint32_t m=-1;
  b[0]=b[0] ^ m;       // b=1....1 if sp>=0. b=0 if sp<0.

  FullRefreshBool(b,c,k,n);

  uint32_t z3[n];
  SecAnd(sp,c,z3,k,n); // if sp>=0, then b=1....1 and we select sp

  for(int i=0;i<n;i++)
    z[i]=z2[i] ^ z3[i];
}
  


void Expand(uint32_t *x,uint32_t *xp,int k,int n2,int n)
{
  for(int i=0;i<n/2;i++)
  {
    uint32_t r=genrand(k);
    xp[2*i]=x[i] ^ r;
    xp[2*i+1]=r;
  }
  if ((n & 1)==1) 
  {
    if (n2==n/2)
      xp[n-1]=0;
    else
      xp[n-1]=x[n2-1];
  }
}

// Goubin's first order conversion from arithmetic to Boolean
// Returns x such that A+r=x xor r
uint32_t GoubinAB(uint32_t A,uint32_t r,int k)
{
  uint32_t G=rand32();
  uint32_t T=G << 1;
  uint32_t x=G ^ r;
  uint32_t O=G & x;
  x=T ^ A;
  G=G ^ x;
  G=G & r;
  O=O ^ G;
  G=T & A;
  O=O ^ G;
  for(int i=1;i<k;i++)
  {
    G=T & r;
    G=G ^ O;
    T=T & A;
    G=G ^ T;
    T=G << 1;
  }
  x=x ^ T;
  return x;
}

    
void ConvertAB(uint32_t *A,uint32_t *z,int k,int n)
{
  if(n==1)
  {
    z[0]=A[0];
    return;
  }

  if(n==2)
  {
    z[0]=GoubinAB(A[0],A[1],k);
    z[1]=A[1];
    return;
  }

  uint32_t x[n/2];
  ConvertAB(A,x,k,n/2);
  uint32_t xp[n];
  Expand(x,xp,k,n/2,n);
  
  uint32_t y[(n+1)/2];
  ConvertAB(A+n/2,y,k,(n+1)/2);
  uint32_t yp[n];
  Expand(y,yp,k,(n+1)/2,n);

  SecAdd(xp,yp,z,k,n);
}



void ConvertABModp(uint32_t *A,uint32_t *z,uint32_t p,int k,int n)
{
  if(n==1)
  {
    z[0]=A[0];
    return;
  }

  uint32_t x[n/2];
  ConvertABModp(A,x,p,k,n/2);
  uint32_t xp[n];
  Expand(x,xp,k,n/2,n);
  
  uint32_t y[(n+1)/2];
  ConvertABModp(A+n/2,y,p,k,(n+1)/2);
  uint32_t yp[n];
  Expand(y,yp,k,(n+1)/2,n);

  SecAddModp(xp,yp,z,p,k,n);
}





   


