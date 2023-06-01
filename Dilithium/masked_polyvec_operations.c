#include "masked_polyvec_operations.h"
#include "masking_interface.h"
#include "../Masking/random.h"



void masked_polyvec_matrix_pointwise_montgomery(masked_polyveck *t, const polyvecl mat[K], const masked_polyvecl *v){
  for(int i=0; i < N_SHARES; ++i) polyvec_matrix_pointwise_montgomery(&(t->shares[i]), mat, &(v->shares[i]));
}

void masked_polyveck_reduce(masked_polyveck *v){
  for(int i=0; i < N_SHARES; ++i) polyveck_reduce(&(v->shares[i])); 
}
void masked_polyveck_invntt_tomont(masked_polyveck *v){
  for(int i=0; i < N_SHARES; ++i) polyveck_invntt_tomont(&(v->shares[i])); 
}

void masked_polyvecl_ntt(masked_polyvecl *v) {
  for(int i=0; i < N_SHARES; ++i) polyvecl_ntt(&(v->shares[i])); 
}

void masked_polyveck_ntt(masked_polyveck *v) {
  for(int i=0; i < N_SHARES; ++i) polyveck_ntt(&(v->shares[i])); 
}

void masked_polyvecl_invntt_tomont(masked_polyvecl *v) {
  for(int i=0; i < N_SHARES; ++i) polyvecl_invntt_tomont(&(v->shares[i]));
}

void masked_polyvecl_pointwise_poly_montgomery(masked_polyvecl *r, const poly *a, const masked_polyvecl *v) {
  for(int i=0; i < N_SHARES; ++i) polyvecl_pointwise_poly_montgomery(&(r->shares[i]), a, &(v->shares[i]));
}

void masked_polyvecl_add(masked_polyvecl *w, const masked_polyvecl *u, const masked_polyvecl *v) {
  for(int i=0; i < N_SHARES; ++i) polyvecl_add(&(w->shares[i]), &(u->shares[i]), &(v->shares[i]));
}

void masked_polyvecl_reduce(masked_polyvecl *v){
  for(int i=0; i < N_SHARES; ++i) polyvecl_reduce(&(v->shares[i]));
}


void masked_polyveck_pointwise_poly_montgomery(masked_polyveck *r, const poly *a, const masked_polyveck *v){
  for(int i=0; i < N_SHARES; ++i) polyveck_pointwise_poly_montgomery(&(r->shares[i]), a, &(v->shares[i]));
}
void masked_polyveck_sub(masked_polyveck *w, const masked_polyveck *u, const masked_polyveck *v){
  for(int i=0; i < N_SHARES; ++i) polyveck_sub(&(w->shares[i]), &(u->shares[i]), &(v->shares[i]));
}

void masked_polyveck_caddq(masked_polyveck *v){
  for(int i=0; i < N_SHARES; ++i) polyveck_caddq(&(v->shares[i]));
}










