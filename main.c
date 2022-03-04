//Comentario
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define L 100
#define T0 20.0
#define TL 100.0
#define TR 100.0
#define ERR 0.001
#define N 1000
#define C 1e-5
#define DELTAT 1
#define DELTAX 0.02

void print(int n, double arrayToPrint[]) {
	for (int i = 0; i < n; i++){
		printf("%f\n", arrayToPrint[i]);
	}
	return;
}

double newTemp(double previousTj, double currentTj, double nextTj){
	double newValue = currentTj + (C*DELTAT/(DELTAX*DELTAX))*(previousTj - 2*currentTj + nextTj);
	//printf("%f - %f - %f\n", previousTj, currentTj, nextTj);
	return newValue;
}

double checkPrecisionRequired(int n, double array[]){
	double diff = 0;
	for(int i = 1; i<n-1; i++){
		diff = diff + abs(array[i] - array[i - 1]);
	}
	return diff/n;
}

int main(int argc, char* argv[]) {
	double err = ERR;
	int n = N;
	if(argc > 5) {
		err = strtol(argv[1], NULL, 10);
		n = strtol(argv[2], NULL, 10);
	}
	
	double answer[n];
	double temp[n];

	if(argc > 5) {
		answer[0] = strtol(argv[4], NULL, 10);
		answer[N-1] = strtol(argv[5], NULL, 10);
		temp[0] = strtol(argv[4], NULL, 10);
		temp[N-1] = strtol(argv[5], NULL, 10);
		for (int temp1 = 1; temp1<N-1; temp1++){
			answer[temp1] = strtol(argv[3], NULL, 10);
		}
  	} else {
		answer[0] = TL;
		answer[N-1] = TR;
		temp[0] = TL;
		temp[N-1] = TR;
		for (int temp1 = 1; temp1<N-1; temp1++){
			answer[temp1] = T0;
		}
	}	

	while(checkPrecisionRequired(n, answer) > err) {
		for(int j = 1; j<n-1; j++){
			temp[j] = newTemp(answer[j-1], answer[j], answer[j+1]);
		}
		for (int k = 1; k < n-1; k++){
			answer[k] = temp[k];
		}
		//printf("New values\n");
		//print(answer);
	}
	printf("Final values\n");
	print(n, answer);
}
