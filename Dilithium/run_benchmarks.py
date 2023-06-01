from os import system, popen
import sys
PATH = "./bench_res.txt"
open(PATH, 'w').close() # clear file
cur = ""
res = ""

VERBOSE_COMPILE = True
REDIRECT = ""
MAX_ORDER = 1


if not VERBOSE_COMPILE:
	REDIRECT=" >/dev/null"
with open(PATH,'a') as f:
  for mode in [2,3,5]:
    f.write("\nMode: "+str(mode))
    for i in range(1, MAX_ORDER+1):
      print("Compiling for masking of order", i,"and mode", mode)
      #system("make clean > /dev/null && make bench ORDER="+str(i)+" RNG="+str(rng)+REDIRECT)
      system("make clean > /dev/null && make bench ORDER="+str(i)+" MODE="+ str(mode)+REDIRECT)
      print("Running tests...", end ='')
      sys.stdout.flush()
      cur = popen("./bench").read()
      print("Writing to "+PATH+" ...")
      f.write(cur)
      res += cur
      print(cur)
      print("Done.")

