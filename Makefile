all:
	g++-5 -O3 -Wall -Itrng-4.19 -Ltrng-4.19/src/.libs HLM_KMP_Flux_length_Test1_Hwang.cpp -o exe_flux1 -std=c++17 -ltrng4 -fopenmp
run:
	./exe_flux1
