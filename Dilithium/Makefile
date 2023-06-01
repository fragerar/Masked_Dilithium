CC = gcc
#CFLAGS += -Wall -Wextra -Wpedantic -Wredundant-decls \
#  -Wshadow -Wpointer-arith -O3 -march=native -mtune=native 
CFLAGS += -O3 -march=native
NISTFLAGS += -Wno-unused-result -O3
SOURCES = sign.c packing.c polyvec.c poly.c ntt.c reduce.c rounding.c
MASKING_LIB_SOURCES = ../Masking/utils.c ../Masking/random.c ../Masking/dilithium_gadgets.c ../Masking/convtable.c ../Masking/convba_2014.c ../Masking/impconvba.c
MASKED_SOURCES = $(MASKING_LIB_SOURCES) masked_sign.c masking_interface.c masked_polyvec_operations.c bench_masked_sign.c
MASKED_HEADERS = masking_interface.h
HEADERS = config.h params.h api.h sign.h packing.h polyvec.h poly.h ntt.h \
  reduce.h rounding.h symmetric.h randombytes.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h
AES_SOURCES = $(SOURCES) fips202.c aes256ctr.c symmetric-aes.c
AES_HEADERS = $(HEADERS) fips202.h aes256ctr.h

ORDER=1
MODE=3
MACRO = -D MASKING_ORDER=$(ORDER)

.PHONY: all speed shared clean

all: simple_tests



simple_tests: simple_test.c randombytes.c $(MASKED_SOURCES) $(MASKED_HEADERS) $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=$(MODE) \
	  -o main $< randombytes.c $(KECCAK_SOURCES) $(MASKED_SOURCES)


clean:

	rm -f main
