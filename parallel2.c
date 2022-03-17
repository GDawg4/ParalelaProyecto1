/* 
 * Universidad del Valle de Guatemala
 * Computación Paralela y Distribuida
 * Proyecto 1
 * Rodrigo Garoz 18102
 * Jose Miguel Castañeda
 * Douglas de León Molina
 * 
 * Argumentos:
 * 1. err - precisión requerida
 * 2. N - número de intervalos
 * 3. T0 - temperatura inicial de la barra
 * 4. TL - temperatura en la frontera izquierda (x=0)
 * 5. TR - temperatura en la frontera derecha (x=L)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <omp.h>

#define DEFAULT_ERR 0.0001
#define DEFAULT_N 10000
#define DEFAULT_T0 20.0
#define DEFAULT_TL 100.0
#define DEFAULT_TR 40.0

#define L 100
#define c 1e-5
#define C 0.5
#define DELTAT 1
#define DELTAX 0.02
#define TESTING 1

int thread_count = 8;

double err = DEFAULT_ERR;
int N = DEFAULT_N;
double T0 = DEFAULT_T0;
double TL = DEFAULT_TL;
double TR = DEFAULT_TR;
double dx;
double dt;
double used_factor;

void print(int n, double arrayToPrint[]);

double newTemp(double previousTj, double currentTj, double nextTj);

double check_convergence(int n, double temp[], double answer[]);


int main(int argc, char* argv[]) {
	double start_t;
	double end_t;

	if(argc > 5) {
		err = strtol(argv[1], NULL, 10);
		N = strtol(argv[2], NULL, 10);
		T0 = strtol(argv[3], NULL, 10);
		TL = strtol(argv[4], NULL, 10);
		TR = strtol(argv[5], NULL, 10);
	} else {
		if (!TESTING) {
			// Input error
			printf("Input error: ");
			scanf("%lf", &err);
			// Input N
			printf("Input N: ");
			scanf("%d", &N);
			// Input T0
			printf("Input T0: ");
			scanf("%lf", &T0);
			// Input TL
			printf("Input TL: ");
			scanf("%lf", &TL);
			// Input TR
			printf("Input TR: ");
			scanf("%lf", &TR);
		}
	}
	start_t = omp_get_wtime();

	dx = (double)L / (double)N;
	dt = (C * dx * dx) / c;
	used_factor = (c*dt)/(dx*dx);
	printf("err: %lf\n", err);
	printf("N: %d\n", N);
	printf("T0: %lf\n", T0);
	printf("TL: %lf\n", TL);
	printf("TR: %lf\n", TR);
	printf("dx: %lf\n", dx);
	printf("dt: %lf\n", dt);

	double answer[N];
	double temp[N];
	// initialize array values
	answer[0] = TL;
	answer[N-1] = TR;
	temp[0] = TL;
	temp[N-1] = TR;
	int i1;
	#pragma omp parallel for private(i1) schedule(static) num_threads(thread_count)
	for (i1=1; i1<N-1; ++i1){
		answer[i1] = T0;
	}

	int converges = 0;
	int j1 = 0;
	while(!converges) {
		#pragma omp parallel for num_threads(thread_count) shared(temp, answer)
		for(j1 = 1; j1<N-1; ++j1){
			temp[j1] = newTemp(answer[j1-1], answer[j1], answer[j1+1]);
		}

		converges = check_convergence(N, temp, answer);

		memcpy( answer, temp, sizeof(answer) );
	}
	// printf("Final values\n");
	// print(N, answer);
	printf("%lf, %lf, %lf, ..., ", answer[0], answer[1], answer[2]);
    	printf("%lf, %lf, %lf, ..., ", answer[N/2], answer[N/2 +1], answer[N/2 +2]);
    	printf("%lf, %lf, %lf\n", answer[N-3], answer[N-2], answer[N-1]);

	end_t = omp_get_wtime();

	printf("Elapsed time: %f\n", (end_t - start_t));

	return 0;
}


void print(int n, double arrayToPrint[]) {
	for (int i = 0; i < n; ++i){
		printf("%lf ", arrayToPrint[i]);
	}
	printf("\n");
	return;
}


double newTemp(double previousTj, double currentTj, double nextTj){
	return currentTj + used_factor*(previousTj - 2*currentTj + nextTj);
	//printf("%f - %f - %f\n", previousTj, currentTj, nextTj);
	//return newValue;
}
double max(double A[], double B[], int i, int j)
{
    double max_val = A[0]-B[0];

    #pragma omp parallel for schedule (static) num_threads(thread_count) shared(A, B) reduction(max:max_val)
    for (int idx = i; idx < j; ++idx)
       max_val = max_val > fabs(A[idx]-B[idx]) ? max_val : fabs(A[idx]-B[idx]);

    return max_val;
}

double check_convergence(int n, double temp[], double answer[]){
	double max_diff = max(temp, answer, 1, N-1);

	if (max_diff < err) {
		printf("max diff: %lf\n", max_diff);
		return 1;
	}
	else
		return 0;
}
