//Comentario
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define L 100
#define T0 20.0
#define TL 100.0
#define TR 100.0
#define ERR 0.01
#define N 1000
#define C 1e-5
#define DELTAT 1
#define DELTAX 0.02

void print(double arrayToPrint[]) {
	for (int i = 0; i < N; i++){
		printf("%f\n", arrayToPrint[i]);
	}
	return;
}

double newTemp(double previousTj, double currentTj, double nextTj){
	double newValue = currentTj + (C*DELTAT/(DELTAX*DELTAX))*(previousTj - 2*currentTj + nextTj);
	//printf("%f - %f - %f\n", previousTj, currentTj, nextTj);
	return newValue;
}

int main(int argc, char* argv[]) {
	double answer[N];
	double temp[N];
	answer[0] = TL;
	answer[N-1] = TR;
	temp[0] = TL;
	temp[N-1] = TR;
	for (int temp1 = 1; temp1<N-1; temp1++){
		answer[temp1] = T0;
	}
	for(int i = 0; i<1000000; i++){
		for(int j = 1; j<N-1; j++){
			temp[j] = newTemp(answer[j-1], answer[j], answer[j+1]);
		}
		for (int k = 1; k < N-1; k++){
			answer[k] = temp[k];
		}
		//printf("New values\n");
		//print(answer);
	}
	printf("Final values\n");
	print(answer);
}
