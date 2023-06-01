// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.

// This is an implementation of some of the algorithms in [CGMZ22]:
// [CGMZ22]  Jean-Sébastien Coron, François Gérard, Simon Montoya, Rina Zeitoun. High-order
//           Table-based Conversion Algorithms and Masking Lattice-based Encryption.
//           IACR Trans. Cryptogr. Hardw. Embed. Syst. 2022(2): 1-40 (2022)


// We implement:
// - Generic Boolean to arithmetic, with l bits as input, and k bits as output. Algorithm 3.
// - Generic arithmetic to Boolean, with l bits as input, and l bits as output. Algorithm 8.
// - Register-based arithmetic to 1-bit Boolean, with 5,6,7 bits as input. Algorithm 12.
// - Optimized Boolean to arithmetic, with l bits as input and k bits as output. We asssume that l is even. Algorithm 5.
// - Arithmetic shift1, first technique with complexity O(2^l*n^3). Algorithm 6.
// - Arithmetic shift2, second technique with complexity O(2^(2l)*n^2). Algorithm 7.
// - Arithmetic to Boolean conversion with k bits as input, by blocks of ell bits. Algorithm 10. Based on shift1 and shift2.
// - Arithmetic to Boolean conversion with k bits as input, bit-by-bit. Algorithm 11. Based on shift2.

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "random.h"
#include "utils.h"






// Boolean to arithmetic, 1 bit as input, and k bits as output.
// This is Algorithm 2 in signatures.pdf, and corresponds to [SPOG19]
void bool2ArithSPOG(uint32_t *x,uint32_t *y,int k,int n)
{
  y[0]=x[0];
  for(int i=1;i<n;i++)
  {
    y[i]=0;
    refreshArith(y,k,i+1);
    
    for(int j=0;j<(i+1);j++)
      y[j]=y[j]*(1-(x[i]<<1));
    y[0]=y[0]+x[i];
  }
  refreshArith(y,k,n);
}



// Boolean to arithmetic, 1 bit as input, and modulo q as ouptput
// This is Algorithm 1 in signatures.pdf, and corresponds to [SPOG19]
// We consider the NI version, without the last LinearRefresh
void bool2ArithSPOGmodq(uint32_t *x,uint32_t *y,int q,int n)
{
  y[0]=x[0];
  for(int i=1;i<n;i++)
  {
    y[i]=0;
    refreshArithModp(y,q,i+1);
    for(int j=0;j<(i+1);j++)
      y[j]=(q+y[j]*(1-2*x[i])) % q;
    y[0]=(y[0]+x[i]) % q;
  }
}


void bool2ArithSPOGmodq64(uint64_t *x, uint64_t *y, uint64_t q, int n)
{
  y[0]=x[0];
  for(int i=1;i<n;i++)
  {
    y[i]=0;
    refreshArithModp64(y,q,i+1);
    for(int j=0;j<(i+1);j++)
      y[j]=(q+y[j]*(1-2*x[i])) % q;
    y[0]=(y[0]+x[i]) % q;
  }
}



// Boolean to arithmetic, l bits as input, and modulo q as ouptput
// This is Algorithm 1 in signatures.pdf, and corresponds to [SPOG19]
// We consider the NI version, without the last LinearRefresh
void bool2ArithSPOGmodqMulti(uint32_t *x,uint32_t *y,int l,int q,int n)
{
  for(int i=0;i<n;i++) y[i]=0;
  for(int j=0;j<l;j++)
  {
    uint32_t z[n];
    for(int i=0;i<n;i++)
      z[i]=(x[i] >> j) & 1;
    uint32_t v[n];
    bool2ArithSPOGmodq(z,v,q,n);
    for(int i=0;i<n;i++)
      y[i]+=(v[i] << j) % q;
  }
}


void FullRefreshArith(uint32_t *a,uint32_t *b,int n)
{
  for(int i=0;i<n;i++)
    b[i]=a[i];

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32();
      b[i]+=tmp;
      b[j]-=tmp;
    }
  }
}


// Arithmetic shift, 1-bit shift with complexity O(n^2)
// This corresponds to Algorithm3 in signatures.pdf
void shift3(uint32_t *x,uint32_t *a,int k,int n)
{
  uint32_t b[n];
  for(int i=0;i<n;i++)
    b[i]=x[i] & 1;

  uint32_t y[n];
  bool2ArithSPOG(b,y,k,n);

  uint32_t u[n];
  for(int i=0;i<n;i++)
    u[i]=x[i]-y[i];

  for(int i=0;i<n;i++)
  {
    u[n-1]+=(u[i] & 1);
    u[i]-=(u[i] & 1);
  }
  for(int i=0;i<n;i++)
    a[i]=u[i] >> 1;
}


 

// Optimized algorithm, doing 1 shift at a time
// This corresponds to [Algorithm11,CGMZ22]
// Based on shift3
void arith2BoolOpti3NI(uint32_t *x,uint32_t *y,int k,int n)
{
  uint32_t a[n];
  for(int i=0;i<n;i++)
    a[i]=x[i];
  
  for(int i=0;i<n;i++)
    y[i]=0;
  
  for(int j=0;j<k;j++)
  {
    for(int i=0;i<n;i++)
      y[i]+=((a[i] & 1) << j);
    uint32_t a2[n];
    if(j<k-1)
    {
      shift3(a,a2,k-j,n);
      for(int i=0;i<n;i++) a[i]=a2[i];
    }
  }
}



