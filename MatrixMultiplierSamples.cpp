#include <iostream>
#include <cstdlib> 
#include <string>
#include <ctime>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

/*hold calculated sample sizes CI 95% MoE 5% for mini samples*/
int sampleSizesSerial[10];
int sampleSizesParallel[10];
int sampleSizesParallelOptimized[10];

//to calculate sample size, average and sd is needed
const int sampleSizeEst = 10;

double serialTimeStamp[sampleSizeEst];
double parallelTimeStamp[sampleSizeEst];
double parallelOptTimeStamp[sampleSizeEst];



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
/*function initialization*/

int main(){

	initializeMatrix(1000);	//initialize matrices with random values
	cout << "Matrix A and Matrix B was initialized: \n";

	for(int i = 0; i<10;i++){


		double averageSer = 0.0;
		double sdSer = 0.0;

		double averagePar = 0.0;
		double sdPar = 0.0;

		for(int j = 0; j<sampleSizeEst;j++){

			clearArray(matrixC,(i+1)*100);
			clearArray(matrixD,(i+1)*100);
			clearArray(matrixE,(i+1)*100);
			serialTimeStamp[j] =  multiplyMatrixSerial((i+1)*100);
			averageSer += serialTimeStamp[j];
			parallelTimeStamp[j] =  multiplyMatrixParallel((i+1)*100);
			averagePar += parallelTimeStamp[j];
		}

		averageSer /= sampleSizeEst;
		averagePar /= sampleSizeEst;

		for(int j = 0; j<sampleSizeEst;j++){
			sdSer += pow((serialTimeStamp[j] - averageSer), 2);
			sdPar += + pow((parallelTimeStamp[j] - averagePar), 2);
		}

		sdSer = sqrt(sdSer);
		sdPar = sqrt(sdPar);

		sampleSizesSerial[i] = (int) pow(((100.0 * 1.96 * (double) sdSer) / (5 * averageSer)),2);
		sampleSizesParallel[i] = (int) pow(((100.0 * 1.96 * (double) sdPar) / (5 * averagePar)),2);

		sampleSizesSerial[i] = max(sampleSizesSerial[i], sampleSizeEst);
		sampleSizesParallel[i] = max(sampleSizesParallel[i], sampleSizeEst);

		cout << "Sample size for mulitply in Serial "<<(i+1)*100<<" = " << sampleSizesSerial[i] <<endl;
		cout << "Sample size for mulitply in Parallel "<<(i+1)*100<<" = " << sampleSizesParallel[i] <<endl<<endl;



	}

	cout << "Calculating sample size is finished" <<endl;

	for(int i = 0; i<10;i++){

		clearArray(matrixC,(i+1)*100);
		clearArray(matrixD,(i+1)*100);
		clearArray(matrixE,(i+1)*100);

		double averageSer = 0.0;

		double averagePar = 0.0;

		for(int j = 0; j<sampleSizesSerial[i];j++){
			clearArray(matrixC,(i+1)*100);
			averageSer +=  multiplyMatrixSerial((i+1)*100);
		}

		for(int j = 0; j<sampleSizesParallel[i];j++){
			clearArray(matrixD,(i+1)*100);
			averagePar +=  multiplyMatrixParallel((i+1)*100);
		}

		averageSer /= sampleSizesSerial[i];
		averagePar /= sampleSizesParallel[i];

		cout << "Average time in sec to mulitply in Serial for "<<(i+1)*100<<" inputs = " << averageSer <<endl;
		cout << "Average time in sec to mulitply in Parallel for "<<(i+1)*100<<" inputs = " << averagePar <<endl<<endl;
	}

	
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

double multiplyMatrixParallel(int n){
	//matrix multiplication
	double begin = omp_get_wtime();	//start time
	#pragma omp parallel for schedule(dynamic)		//this for loop will be evenly distributed among n number of threads threads and be calculated independent
	for (int i = 0; i < n; i++){
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