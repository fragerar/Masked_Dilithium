#include "masking_interface.h"
#include "../Masking/dilithium_gadgets.h"
#include "params.h"
#include "polyvec.h"



int32_t canon_to_centered(uint32_t x){
  int32_t res = (int32_t) x;
  if (res >= Q/2) res -= Q;
  return res;
}

int32_t center(int32_t x){
  /* Maps elements of Z_q in [-Q, ..., Q] to representatives in [-Q/2, ..., Q/2[ */
  x += Q;
  x %= Q;
  if (x > Q/2) x -= Q;
  return x;
}

uint32_t center_to_canon(int32_t x){
  /* Maps elements of Z_q in [-Q/2, ..., Q/2[ to representatives in  [0, ..., Q[*/
  if (x < 0) x += Q; 
  return x;
}



void masked_sample_y(masked_polyvecl* masked_y){
  polyvecl* y;
  uint32_t masked_coeff[N_SHARES];

  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      gen_y(masked_coeff, N_SHARES);
      for(int k=0; k < N_SHARES; ++k){
          y = &(masked_y->shares[k]);
          y->vec[i].coeffs[j] = canon_to_centered(masked_coeff[k]);
      }
    }
  }
}




int masked_rejection_sampling_z(masked_polyvecl* mz){
  /* Should return 1 if reject */
  uint32_t temp[N_SHARES];

  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      for(int k=0; k < N_SHARES; ++k){
        temp[k] = center_to_canon(mz->shares[k].vec[i].coeffs[j]);
      }
    if (ABC_rejection_sampling(temp, 0, N_SHARES) == 0) return 1;
    }
  }
  return 0;
}

int masked_rejection_sampling_r(masked_polyveck* mr){
  /* Should return 1 if reject */
  uint32_t temp[N_SHARES];

  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      for(int k=0; k < N_SHARES; ++k){
        temp[k] = center_to_canon(mr->shares[k].vec[i].coeffs[j]);
      }
    if (ABC_rejection_sampling(temp, 1, N_SHARES) == 0) return 1;
    }
  }
  return 0;
}



void masked_decompose(polyveck* r1, masked_polyveck* mr0, masked_polyveck* mr){


  uint32_t temp_r[N_SHARES], temp_r0[N_SHARES];


  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      for(int k=0; k < N_SHARES; ++k) temp_r[k] = mr->shares[k].vec[i].coeffs[j];
      decompose2(temp_r, &(r1->vec[i].coeffs[j]), temp_r0, N_SHARES);
      for(int k=0; k < N_SHARES; ++k) mr0->shares[k].vec[i].coeffs[j] = temp_r0[k];
    }
  } 


}




void unmask_polyvecl(masked_polyvecl* mpv, polyvecl* pv){
  int32_t temp;
  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      temp = 0;
      for(int k=0; k < N_SHARES; ++k) temp = (temp + (mpv->shares[k]).vec[i].coeffs[j])%Q;
      pv->vec[i].coeffs[j] = canon_to_centered((temp+2*Q)%Q);
    }
  }
}


void unmask_polyveck(masked_polyveck* mpv, polyveck* pv){
  int32_t temp;
  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      temp = 0;
      for(int k=0; k < N_SHARES; ++k) temp = (temp + (mpv->shares[k]).vec[i].coeffs[j])%Q;
      pv->vec[i].coeffs[j] = canon_to_centered((temp+2*Q)%Q);
    }
  }
}



void mask_polyvecl(masked_polyvecl* mpv, polyvecl* pv){
 /* Takes an unmasked polyvecl and mask it (arith mod q), should not be used outside of testing*/
  int32_t temp;
  for(int i=0; i < L; ++i){
    for(int j=0; j < N; ++j){
      ((mpv->shares[0]).vec[i]).coeffs[j] = (pv->vec[i]).coeffs[j];
    }
  }

  for(int k=1; k < N_SHARES; ++k){
    for(int i=0; i < L; ++i){
      for(int j=0; j < N; ++j){
        temp = ((int32_t)rand32())%Q;
        ((mpv->shares[k]).vec[i]).coeffs[j] = center(temp);
        ((mpv->shares[0]).vec[i]).coeffs[j] = center((((mpv->shares[0]).vec[i]).coeffs[j] - temp)%Q);
      }
    }
  }
}


void mask_polyveck(masked_polyveck* mpv, polyveck* pv){
 /* Takes an unmasked polyveck and mask it (arith mod q), should not be used outside of testing*/
  int32_t temp;
  for(int i=0; i < K; ++i){
    for(int j=0; j < N; ++j){
      ((mpv->shares[0]).vec[i]).coeffs[j] = (pv->vec[i]).coeffs[j];
    }
  }

  for(int k=1; k < N_SHARES; ++k){
    for(int i=0; i < K; ++i){
      for(int j=0; j < N; ++j){
        temp = ((int32_t)rand32())%Q;
        ((mpv->shares[k]).vec[i]).coeffs[j] = center(temp);
        ((mpv->shares[0]).vec[i]).coeffs[j] = center((((mpv->shares[0]).vec[i]).coeffs[j] - temp)%Q);
      }
    }
  }
}

void print_masked_polyvecl(masked_polyvecl* mpv){
  polyvecl t;
  unmask_polyvecl(mpv, &t);
  for(int i=0; i < L; ++i){
    printf("Pv[%i]: ", i);
    for(int j=0; j < 8; ++j) printf("%i ", t.vec[i].coeffs[j]);
    printf("\n");
  }
  printf("\n");
}

void print_polyvecl(polyvecl* t){
  for(int i=0; i < L; ++i){
    printf("Pv[%i]: ", i);
    for(int j=0; j < 8; ++j) printf("%i ", t->vec[i].coeffs[j]);
    printf("\n");
  }
  printf("\n");
}
