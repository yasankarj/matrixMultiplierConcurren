#include <iostream>
#include <cstdlib> 
#include <string>
#include <ctime>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

//to calculate sample size, average and sd is needed
const int noAttempts = 10;		//for population of 100, this is the sample size for CI = 95% and ME = 5%

/*matrix declaration*/
double matrixA[1000][1000];
double matrixB[1000][1000];
double matrixC[1000][1000];
double matrixD[1000][1000];
double matrixE[1000][1000];
/*matrix declaration*/

/*function initialization*/
void initializeMatrix(int n);
void clearArray(double matrix[1000][1000], int n);
void writeMatrix(double matrix[1000][1000], char const* fileName ,int n);
double multiplyMatrixSerial(int n);
double multiplyMatrixParallel(int n);
double multiplyMatrixParallelOpt(int n);
/*function initialization*/

int main(){

	int N = 0;

	cout <<"Enter dimensions"<<endl;
	cin >> N;

	initializeMatrix(1000);	//initialize matrices with random values
	cout << "Matrix A and Matrix B was initialized: \n";

	double averageSer = 0.0;
	double averagePar = 0.0;
	double averageParOpt = 0.0;

	for(int j = 0; j<noAttempts;j++){
		clearArray(matrixC,N);
		clearArray(matrixD,N);
		clearArray(matrixE,N);
		averageSer +=  multiplyMatrixSerial(N);
		averagePar +=  multiplyMatrixParallel(N);
		averageParOpt += multiplyMatrixParallelOpt(N);
	}

	averageSer /= noAttempts;
	averagePar /= noAttempts;
	averageParOpt /= noAttempts;

	cout << "Average time in sec to mulitply in Serial for "<<N<<" inputs = " << averageSer <<endl;
	cout << "Average time in sec to mulitply in Parallel for "<<N<<" inputs = " << averagePar <<endl;
	cout << "Average time in sec to mulitply in Parallel optimized for "<<N<<" inputs = " << averageParOpt << endl <<endl;
	
}

void initializeMatrix(int n){
	//initialize both matrices with random values
	for (int i=0; i<n; i++){		
		for (int j=0; j<n; j++){
			matrixA[i][j] = rand();
			matrixB[i][j] = rand();
		}
	}
}

double multiplyMatrixParallelOpt(int n){
	//matrix multiplication
	double begin = omp_get_wtime();	//start time
	#pragma omp parallel for
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				matrixE[i][j] += matrixA[i][k] * matrixB[k][j];
			}
		}
	}
	double end = omp_get_wtime();;		//end time

	return (double)(end-begin);		//duration for multiplication is calculated
}


double multiplyMatrixParallel(int n){
	//matrix multiplication
	double begin = omp_get_wtime();	//start time
	#pragma omp parallel for schedule(dynamic)		//this for loop will be evenly distributed among n number of threads threads and be calculated independent
	for (int i = 0; i < n; i++){
		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				matrixD[i][j] += matrixA[i][k] * matrixB[k][j];
			}
		}
	}
	double end = omp_get_wtime();;		//end time

	return (double)(end-begin);		//duration for multiplication is calculated
}

double multiplyMatrixSerial(int n){
	//matrix multiplication
	double begin = omp_get_wtime();	//start time
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
			}
		}
	}
	double end = omp_get_wtime();	//end time

	return (double)(end-begin);		//duration for multiplication is calculated
}

void clearArray(double matrix[1000][1000], int n){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				matrix[i][j] = 0.0;
			}
		}
	}
}