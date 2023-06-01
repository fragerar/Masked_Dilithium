#ifndef RANDOM_H
#define RANDOM_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#ifndef RNG_MODE
#define RNG_MODE 1
#endif

#define COUNT
#ifdef COUNT
extern uint64_t count_rand;
#endif


uint16_t rand16();
uint32_t rand32();
uint64_t rand64();
#endif
