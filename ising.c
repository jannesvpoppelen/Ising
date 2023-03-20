#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "xoroshiro.h"

#define N 50
#define SIZE (N*N)
#define J 1
#define Tmin 1.9
#define Tmax 2.5
#define dT 0.01

int **spin;
double avgM=0.0,avgE=0.0;  //Initialize ensemble averages

static inline double rng(){return (double)next() / (double)UINT64_MAX;}

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

static inline char sign(int n) { /* Used in print_config */
    return n >= 0 ? '+' : '-';
} 

static void print_config(int** s){	/* Prints spin configuration. */
	printf("Spin configuration\n");
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			printf("%c%d ",sign(s[i][j]),1);
		}
		printf("\n");
	}
}

void allocate_spin(int*** s,int M) {  /* Allocate memory of 2D array. Inspired by HPP Labs. */
    *s=(int**)malloc(M*sizeof(int*));
    for (int i = 0; i < M; i++) {
        (*s)[i] = (int *)malloc((M) * sizeof(int));
    }
}

void free_spin(int **s, int M) {	/* Free memory of 2D array. Inspired by HPP Labs. */
    for (int i = 0; i < M; i++) {
        free(s[i]);
    }
    free(s);
}

static inline double magnetization(int** s){	/* Calculates magnetization M. */
	double M=0;
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			 M+=s[i][j];
		}
	}
	return (double)(M/SIZE);
}
	
static inline double energy(int** s){	/* Calculates energy E. */
	double E=0;
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			E+= -J*(s[i][j] * (s[i-1][j] + s[i+1][j] + s[i][j-1] + s[i][j+1])); //Hamiltonian
		}
	}
return (double) E/(2*SIZE);
}	

static inline double deltaE(int** s,int i,int j){	/* Calculates change in energy dE as a result of flipped spin. */
		double dE=(double)2*J*s[i][j]*(s[i][j-1]+s[i-1][j]+s[i+1][j]+s[i][j+1]); //Factor of two for contributing through all pairs.
return dE;
}

static void metropolis(double T, int n_thrm, int n_measure,int N_threads){	/* Metropolis sampling to simulate flipping spins. */
	int i,j,thread_id,startcol,endcol;
	double dE,x;
	double beta=1/T; //Loop invariant
	int subsize = N / N_threads;
	int remainder = N%N_threads;
	int thrmsteps= n_thrm/N_threads;
	int measuresteps= n_measure/N_threads;	
	
	#pragma omp parallel num_threads(N_threads) private(thread_id, i, j, dE, startcol, endcol) shared(spin)
	{
	thread_id=omp_get_thread_num();
	startcol=thread_id*subsize+1;
	endcol=startcol+subsize-1+(thread_id==N_threads-1)*remainder;
	thrmsteps=(int) thrmsteps*(1+(thread_id==N_threads-1)*(double)(remainder*N)/SIZE); //divide work among threads
	for(int k=0;k<thrmsteps;k++){
	i=next()%N + 1; //random i between 1 and N
	j=next()% (endcol-startcol+1) + startcol; //random j between startcol and endcol
	dE=deltaE(spin, i,j);
	#pragma atomic
	{
	spin[i][j] *= (dE<=0||exp(-dE*beta)>rng())? -1:1;//-(dE <= 0 || exp(-dE*beta) > x) + (dE > 0 && exp(-dE*beta) <= x)	; //If dE<=0, flip spin, else only flip with Boltzmann factor
	}
	}
	#pragma omp barrier
	}

	#pragma omp parallel num_threads(N_threads) private(thread_id, i, j, dE, startcol, endcol) shared(spin) reduction(+:avgM,avgE)
	{
	thread_id=omp_get_thread_num();
	startcol=thread_id*subsize+1;
	endcol=startcol+subsize-1+(thread_id==N_threads-1)*remainder;
	measuresteps=(int) measuresteps*(1+(thread_id==N_threads-1)*(double)(remainder*N)/SIZE);  //divide work among threads
	for(int k=0;k<measuresteps;k++){
	i=next()%N + 1; //random i between 1 and N
	j=next()%(endcol-startcol+1) + startcol; //random j between startcol and endcol
	dE=deltaE(spin, i,j);
	#pragma atomic
	{
	spin[i][j] *= (dE<=0||exp(-dE*beta)>rng())? -1 : 1; //If dE<=0, flip spin, else only flip with Boltzmann factor
	}
	avgM+=magnetization(spin);	
	avgE+=energy(spin);
	}
	#pragma omp barrier
	}
}

void init(int** s){	/* Sets up individual spins in a random configuration. */
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
		s[i][j] = (rng() <= 0.7) * 2 - 1;	//Cheat the thermalization, set 0.7 to something between 0.6 and 0.7 in case results are unsatisfactory
		}
	}
	//Dealing with periodic boundary conditions.  Stick the array to each other just like sticking a paper onto itself.	
	for (int i = 1; i <= N; i++) {
    s[i][0] = s[i][N];
    s[i][N+1] = s[i][1];
	}
	for (int j = 1; j <= N; j++) {
    s[0][j] = s[N][j];
    s[N+1][j] = s[1][j];
	}
}

int main(int argc, const char *argv[]){
	if(argc!=4){
		printf("Incorrect number of arguments\n");
		printf("Supply Number of thermalization sweeps and number of measure sweeps for the Monte Carlo sampling, Number of threads\n");
	return 0;}
	const int n_thrm=atoi(argv[1])*SIZE;
	const int n_measure=atoi(argv[2])*SIZE;
	int N_threads=atoi(argv[3]);
	printf("N=%d\n",N);
	printf("SIZE=%d\n",SIZE);
	double tijd=get_wall_seconds();
	allocate_spin(&spin, N+2);//Need two extra rows and columns to deal with periodic boundary conditions
	init(spin);
	
	FILE *file; //Write T, E, M
	file=fopen("data.txt","w");
	if (file == NULL) {
        printf("Could not open file\n");
        return 0;
    }
	for(double T=Tmin;T<=Tmax;T+=dT){
	printf("T=%.3f\n",T);
	init(spin);
	metropolis(T,n_thrm,n_measure,N_threads);
	avgE/=(n_measure);
	avgM/=(n_measure);
	fprintf(file, "%.5f %.16f %.16f\n", T, avgE, fabs(avgM));
	}
	fclose(file);
	free_spin(spin,N+2);
	printf("Ising sim took %.4lf wall seconds\n",get_wall_seconds()-tijd);
	return 0;}