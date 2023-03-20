ising: ising.c
	gcc -c xoroshiro.c -Ofast -mavx 
	gcc -c ising.c -Ofast -mavx -fopenmp -g -lm
	gcc ising.o xoroshiro.o -Ofast -mavx -fopenmp -g -lm -o ising
	
clean: ising xoroshiro.o
	rm ising.o
	rm xoroshiro.o
	rm ising