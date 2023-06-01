#ifndef CONVBA2014_H
#define CONVBA2014_H

void ConvertAB(uint32_t *A,uint32_t *z,int k,int n);
void ConvertABModp(uint32_t *A,uint32_t *z,uint32_t p,int k,int n);


void impconvBA(uint32_t *D_,uint32_t *x,int n); //[BCZ18]
void impconvBA64(uint64_t *D_,uint64_t *x,int n);
#endif