/* 
 * Universidad del Valle de Guatemala
 * Computación Paralela y Distribuida
 * Proyecto 1
 * Rodrigo Garoz
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

#define DEFAULT_ERR 0.001
#define DEFAULT_N 10000
#define DEFAULT_T0 20.0
#define DEFAULT_TL 100.0
#define DEFAULT_TR 70.0

#define L 100
#define c 1e-5
#define C 0.5
#define DELTAT 1
#define DELTAX 0.02
#define TESTING 1


double err = DEFAULT_ERR;
int N = DEFAULT_N;
double T0 = DEFAULT_T0;
double TL = DEFAULT_TL;
double TR = DEFAULT_TR;
double dx;
double dt;


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
	for (int i=1; i<N-1; i++){
		answer[i] = T0;
	}

	int converges = 0;
	while(!converges) {
        #pragma omp parallel for shared(answer, temp)
		for(int j = 1; j<N-1; j++){
			temp[j] = newTemp(answer[j-1], answer[j], answer[j+1]);
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
	for (int i = 0; i < n; i++){
		printf("%lf ", arrayToPrint[i]);
	}
	printf("\n");
	return;
}


double newTemp(double previousTj, double currentTj, double nextTj){
	double newValue = currentTj + (c*dt/(dx*dx))*(previousTj - 2*currentTj + nextTj);
	//printf("%f - %f - %f\n", previousTj, currentTj, nextTj);
	return newValue;
}


double check_convergence(int n, double temp[], double answer[]){
	double diff = 0.0;
	double max_diff = 0.0;

    #pragma omp parallel for shared(max_diff) private(diff)
	for(int i = 1; i<N-1; i++){
		diff = fabs(temp[i] - answer[i]);
		if (diff > max_diff) {
            #pragma omp critical
            max_diff = diff;
        }
	}

	if (max_diff < err) {
		printf("max diff: %lf\n", max_diff);
		return 1;
	} 
	else
		return 0;
}
