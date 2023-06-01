// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.

// This is an implementation of the Boolean to arithmetic conversion algorithm from [BCZ18]:
// [BCZ18] Luk Bettale, Jean-Sébastien Coron, Rina Zeitoun. Improved High-Order Conversion
//         From Boolean to Arithmetic Masking. IACR Trans. Cryptogr. Hardw. Embed. Syst. 2018(2): 22-45 (2018)

// We also implement the variant algorithm working for any integer q, not only a power-of-two.

//#define DEBUG

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "random.h"
#include "convtable.h"
#include "utils.h"
#include "convba_2014.h"




uint32_t Psi(uint32_t x,uint32_t y)
{
  return (x ^ y)-y;
}

uint32_t Psi0(uint32_t x,uint32_t y,int n)
{
  return Psi(x,y) ^ ((~n & 1) * x);
}

void copy(uint32_t *x,uint32_t *y,int n)
{
  for(int i=0;i<n;i++) x[i]=y[i];
}

static void impconvBA_rec(uint32_t *D_,uint32_t *x,int n);

void impconvBA(uint32_t *D_,uint32_t *x,int n)
{
  uint32_t x_ext[n+1];
  copy(x_ext,x,n);
  x_ext[n] = 0;
  impconvBA_rec(D_, x_ext, n);
}

// here, x contains n+1 shares
static void impconvBA_rec(uint32_t *D_,uint32_t *x,int n)
{  
  if (n==2)
  {
    uint32_t r1=rand32();
    uint32_t r2=rand32();
    uint32_t y0=(x[0] ^ r1) ^ r2;
    uint32_t y1=x[1] ^ r1;
    uint32_t y2=x[2] ^ r2;
    
    uint32_t z0=y0 ^ Psi(y0,y1);
    uint32_t z1=Psi(y0,y2);
    
    D_[0]=y1 ^ y2;
    D_[1]=z0 ^ z1;

    #ifdef DEBUG
      assert((x[0] ^ x[1] ^ x[2])==(D_[0] + D_[1]));
    #endif
     
    return;
  }

  uint32_t y[n+1];
  copy(y,x,n+1);

  refreshBool(y,32,n+1);

  uint32_t z[n];

  z[0]=Psi0(y[0],y[1],n);
  for(int i=1;i<n;i++)
    z[i]=Psi(y[0],y[i+1]);

  #ifdef DEBUG
  assert(xorop(x,n+1)==(xorop(y+1,n)+xorop(z,n)));
  #endif

  uint32_t A[n-1],B[n-1];
  impconvBA_rec(A,y+1,n-1);
  impconvBA_rec(B,z,n-1);
  
  for(int i=0;i<n-2;i++)
    D_[i]=A[i]+B[i];

  D_[n-2]=A[n-2];
  D_[n-1]=B[n-2];

  #ifdef DEBUG
  assert(xorop(x,n+1)==addop(D_,n));
  #endif
}




uint64_t Psi64(uint64_t x,uint64_t y)
{
  return (x ^ y)-y;
}

uint64_t Psi064(uint64_t x,uint64_t y,int n)
{
  return Psi64(x,y) ^ ((~n & 1) * x);
}

void copy64(uint64_t *x,uint64_t *y,int n)
{
  for(int i=0;i<n;i++) x[i]=y[i];
}

static void impconvBA_rec64(uint64_t *D_,uint64_t *x,int n);

void impconvBA64(uint64_t *D_,uint64_t *x,int n)
{
  uint64_t x_ext[n+1];
  copy64(x_ext,x,n);
  x_ext[n] = 0;
  impconvBA_rec64(D_, x_ext, n);
}

// here, x contains n+1 shares
static void impconvBA_rec64(uint64_t *D_,uint64_t *x,int n)
{  
  if (n==2)
  {
    uint64_t r1=rand64();
    uint64_t r2=rand64();
    uint64_t y0=(x[0] ^ r1) ^ r2;
    uint64_t y1=x[1] ^ r1;
    uint64_t y2=x[2] ^ r2;
    
    uint64_t z0=y0 ^ Psi64(y0,y1);
    uint64_t z1=Psi64(y0,y2);
    
    D_[0]=y1 ^ y2;
    D_[1]=z0 ^ z1;

    #ifdef DEBUG
      assert((x[0] ^ x[1] ^ x[2])==(D_[0] + D_[1]));
    #endif
     
    return;
  }

  uint64_t y[n+1];
  copy64(y,x,n+1);

  refreshBool64(y,n+1);

  uint64_t z[n];

  z[0]=Psi064(y[0],y[1],n);
  for(int i=1;i<n;i++)
    z[i]=Psi64(y[0],y[i+1]);

  #ifdef DEBUG
  assert(xorop(x,n+1)==(xorop(y+1,n)+xorop(z,n)));
  #endif

  uint64_t A[n-1],B[n-1];
  impconvBA_rec64(A,y+1,n-1);
  impconvBA_rec64(B,z,n-1);
  
  for(int i=0;i<n-2;i++)
    D_[i]=A[i]+B[i];

  D_[n-2]=A[n-2];
  D_[n-1]=B[n-2];

  #ifdef DEBUG
  assert(xorop(x,n+1)==addop(D_,n));
  #endif
}



