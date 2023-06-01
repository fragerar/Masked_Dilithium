#ifndef UTILS_H
#define UTILS_H


int64_t cpucycles(void);
uint32_t genrand(int l);

void share(uint32_t x,uint32_t a[],int n);
void share64(uint64_t x,uint64_t a[],int n);
void refreshBool(uint32_t a[],int l,int n);
void refreshBool64(uint64_t a[],int n);
void FullRefreshBool(uint32_t *a,uint32_t *b,int l,int n);
void refreshArith(uint32_t a[],int l,int n);
void refreshArith64(uint64_t a[],int l,int n);
void refreshArithModp(uint32_t a[],uint32_t p,int n);
void refreshArithModp64(uint64_t a[],uint64_t p,int n);

uint32_t xorop(uint32_t a[],int n);
uint32_t addop(uint32_t a[],int l,int n);
uint32_t addopmodp(uint32_t *a,int p,int n);


void SecAnd(uint32_t *a,uint32_t *b,uint32_t *c,int k,int n);
void SecAdd(uint32_t *x,uint32_t *y,uint32_t *z,int k,int n);

void SecMul(uint32_t *a,uint32_t *b,uint32_t *c,int n);
void SecMultModp(uint32_t *a,uint32_t *b,uint32_t *c,uint32_t p,int n);



void printShares(uint32_t *a,int n);
void printArith(uint32_t *a, int p, int n);
void printBool(uint32_t *a, int n);

void printArith64(uint64_t *a, uint64_t p, int n);
void printBool64(uint64_t *a, int n);



#endif