#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"

#include "masked_sign.h"
#include "masking_interface.h"
#include "masked_polyvec_operations.h"


int masked_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen, uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0)
{
  unsigned int n;
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  masked_polyvecl my, mz;
  masked_polyveck mw1, mh, mw0;
  polyvecl mat[K], y, z;
  polyveck w1, w0, h;
  poly cp;
  keccak_state state;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + CRHBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;

  

  /* Compute CRH(tr, msg) */
  shake256_init(&state);
  shake256_absorb(&state, tr, CRHBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

//#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
/* #else
  crh(rhoprime, key, SEEDBYTES + CRHBYTES);
#endif */

  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);


/*   polyvecl_ntt(&s1);
  polyveck_ntt(&s2); */
  masked_polyvecl_ntt(ms1);
  masked_polyveck_ntt(ms2);
  polyveck_ntt(t0);

  int iter = 0;
rej:
  ++iter;
  
  /* Sample intermediate vector y */
  //polyvecl_uniform_gamma1(&y, rhoprime, nonce++);
  masked_sample_y(&my);
  //z = y;
  mz = my;
  //polyvecl_ntt(&z);
  masked_polyvecl_ntt(&mz);

  /* Matrix-vector multiplication */
/*   polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1); */
  masked_polyvec_matrix_pointwise_montgomery(&mw1, mat, &mz);
  masked_polyveck_reduce(&mw1);
  masked_polyveck_invntt_tomont(&mw1);



  /* Decompose w and call the random oracle */
  //polyveck_caddq(&w1);
  //polyveck_decompose(&w1, &w0, &w1);
  masked_polyveck_caddq(&mw1);
  masked_decompose(&w1, &mw0, &mw1);  

  
  polyveck_pack_w1(sig, &w1);

  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);
  poly_ntt(&cp);

  /* Compute z, reject if it reveals secret */
  //polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  //polyvecl_invntt_tomont(&z);
  //polyvecl_add(&z, &z, &y);
  //polyvecl_reduce(&z);

  masked_polyvecl_pointwise_poly_montgomery(&mz, &cp, ms1);  
  masked_polyvecl_invntt_tomont(&mz);
  masked_polyvecl_add(&mz, &mz, &my);
  masked_polyvecl_reduce(&mz);

  if (masked_rejection_sampling_z(&mz))
    goto rej;




  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  //polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  //polyveck_invntt_tomont(&h);
  //polyveck_sub(&w0, &w0, &h);
  //polyveck_reduce(&w0);

  masked_polyveck_pointwise_poly_montgomery(&mh, &cp, ms2);
  masked_polyveck_invntt_tomont(&mh);
  masked_polyveck_sub(&mw0, &mw0, &mh);
  masked_polyveck_reduce(&mw0);
  if (masked_rejection_sampling_r(&mw0))
    goto rej;


/*   if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    goto rej;

  if(polyveck_chknorm(&w0, GAMMA2 - BETA))
    goto rej; */
  

/*   if (masked_rejection_sampling(&mz, &mw0))
    goto rej; */

  unmask_polyvecl(&mz, &z);
  unmask_polyveck(&mw0, &w0);



  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h, &cp, t0);
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
  if(polyveck_chknorm(&h, GAMMA2))
    goto rej;

  polyveck_add(&w0, &w0, &h);
  polyveck_caddq(&w0);
  n = polyveck_make_hint(&h, &w0, &w1);
  if(n > OMEGA)
    goto rej;

  /* Write signature */
  pack_sig(sig, sig, &z, &h);
  *siglen = CRYPTO_BYTES;
  return iter;
}



int masked_crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0)
{
  size_t i;
  int rej;
  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  rej=masked_crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, seedbuf, ms1, ms2, t0);
  *smlen += mlen;
  return rej;
}