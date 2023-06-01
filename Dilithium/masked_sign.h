
#ifndef MASKED_SIGN_H
#define MASKED_SIGN_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include "masking_interface.h"

int masked_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen, uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0);


int masked_crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0);



int bench_masked_crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0, uint64_t* bench_vector);


int bench_masked_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen, uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0, uint64_t* bench_vector);

#endif