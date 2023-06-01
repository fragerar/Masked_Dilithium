COMMON_SOURCES = utils.c random.c
FLAGS = -Wall -O3 -march=native -Wunused-function


dilithium_gadgets: $(COMMON_SOURCES)  dilithium_gadgets.c convtable.c convba_2014.c impconvba.c
	gcc $(FLAGS) -D TESTS_SAMPLE_Y $(COMMON_SOURCES) dilithium_gadgets.c convtable.c convba_2014.c impconvba.c -o main


clean:
	rm -f main 