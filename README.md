# Masked_Dilithium
To replicate Tables 5, 6, 7 and 10: 
1. go to ./Masking/
2. make
3. ./main

To replicate Table 8: 
1. go to ./Dilithium
2. Edit params.h and/or masking_interface.h to chose security level and masking order
3. make
4. ./main
The vector of values output correspond to the cycles for: [NTT, SampleY, AY, Decompose, Z=y+cs1, Reject, w-cs2]  

To automatically run signature benchmarks at multiple orders and for multiple security levels, one can also use the run_benchmarks Python script. 
