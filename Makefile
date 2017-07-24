all:
	g++-5 -O3 -Wall -Itrng-4.19 -Ltrng-4.19/src/.libs HLM_KMP_Flux_Hwang.cpp -o exe_flux -std=c++11 -ltrng4 -fopenmp
run:
	./exe_flux
