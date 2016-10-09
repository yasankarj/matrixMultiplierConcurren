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
int sampleSizesParallelOp[10];

//to calculate sample size, average and sd is needed
const int sampleSizeEst = 10;

double serialTimeStamp[sampleSizeEst];
double parallelTimeStamp[sampleSizeEst];
double parallelTimeStampOp[sampleSizeEst];


/*matrix declaration*/
double matrixA[1000][1000];
double matrixB[1000][1000];
double matrixC_Serial[1000][1000];
double matrixC_Parallel[1000][1000];
double matrixC_ParallelOp[1000][1000];
double matrixD[1000][1000];
double matrixE[1000][1000];
double matrixT[1000][1000];
double matrixTransB[1000][1000];
/*matrix declaration*/

/*function initialization*/
void initializeMatrix(int n);
void clearArray(double matrix[1000][1000], int n);
double multiplyMatrixParallelOp(int nSamples, int matrixSize);
double multiplyMatrixParallel(int matrixSize);
double multiplyMatrixSerial(int matrixSize);
void pMultiply(int n, double matrixA[1000][1000], double matrixB[1000][1000], double matrixC[1000][1000],
               double matrixT[1000][1000]);


int main() {

    initializeMatrix(1000);    //initialize matrices with random values
    cout << "Matrix A and Matrix B was initialized: \n";

    for (int i = 0; i < 10; i++) {


        double averageSer = 0.0;
        double sdSer = 0.0;

        double averagePar = 0.0;
        double sdPar = 0.0;

        double averageParOp = 0.0;
        double sdParOp = 0.0;

        for (int j = 0; j < sampleSizeEst; j++) {

            clearArray(matrixC_Serial, (i + 1) * 100);
            clearArray(matrixC_Parallel, (i + 1) * 100);
            clearArray(matrixC_ParallelOp, (i + 1) * 100);
            clearArray(matrixD, (i + 1) * 100);
            clearArray(matrixE, (i + 1) * 100);
            serialTimeStamp[j] = multiplyMatrixSerial((i + 1) * 100);
            averageSer += serialTimeStamp[j];
            parallelTimeStamp[j] = multiplyMatrixParallel((i + 1) * 100);
            averagePar += parallelTimeStamp[j];
            parallelTimeStampOp[j] = multiplyMatrixParallelOp(sampleSizeEst, (i + 1) * 100);
            averageParOp += parallelTimeStampOp[j];

        }

        averageSer /= sampleSizeEst;
        averagePar /= sampleSizeEst;
        averageParOp /= sampleSizeEst;

        for (int j = 0; j < sampleSizeEst; j++) {
            sdSer += pow((serialTimeStamp[j] - averageSer), 2);
            sdPar += pow((parallelTimeStamp[j] - averagePar), 2);
            sdParOp += pow((parallelTimeStampOp[j] - averageParOp), 2);
        }

        sdSer = sqrt(sdSer);
        sdPar = sqrt(sdPar);
        sdParOp = sqrt(sdParOp);

        sampleSizesSerial[i] = (int) pow(((100.0 * 1.96 * (double) sdSer) / (5 * averageSer)), 2);
        sampleSizesParallel[i] = (int) pow(((100.0 * 1.96 * (double) sdPar) / (5 * averagePar)), 2);
        sampleSizesParallelOp[i] = (int) pow(((100.0 * 1.96 * (double) sdParOp) / (5 * averageParOp)), 2);

        sampleSizesSerial[i] = max(sampleSizesSerial[i], sampleSizeEst);
        sampleSizesParallel[i] = max(sampleSizesParallel[i], sampleSizeEst);
        sampleSizesParallelOp[i] = max(sampleSizesParallelOp[i], sampleSizeEst);

        cout << "Serial" << "SD " << sdSer << " AVG " << averageSer << endl;
        cout << "Parallel" << "SD " << sdPar << " AVG " << averagePar << endl;
        cout << "Optimized Parallel" << "SD " << sdParOp << " AVG " << averageParOp << endl;
        cout << "Sample size for mulitply in Serial " << (i + 1) * 100 << " = " << sampleSizesSerial[i] << endl;
        cout << "Sample size for mulitply in Parallel " << (i + 1) * 100 << " = " << sampleSizesParallel[i] << endl;
        cout << "Sample size for mulitply in Optimized Parallel " << (i + 1) * 100 << " = " <<
        sampleSizesParallelOp[i] << endl << endl;
    }

    cout << "Calculating sample size is finished" << endl;

    for (int i = 0; i < 10; i++) {

        clearArray(matrixC_Serial, (i + 1) * 100);
        clearArray(matrixC_Parallel, (i + 1) * 100);
        clearArray(matrixC_ParallelOp, (i + 1) * 100);

        double averageSer = 0.0;

        double averagePar = 0.0;

        double averageParOp = 0.0;

        for (int j = 0; j < sampleSizesSerial[i]; j++) {
            clearArray(matrixC_Serial, (i + 1) * 100);
            averageSer += multiplyMatrixSerial((i + 1) * 100);
        }

        for (int j = 0; j < sampleSizesParallel[i]; j++) {
            clearArray(matrixC_Parallel, (i + 1) * 100);
            averagePar += multiplyMatrixParallel((i + 1) * 100);
        }

        for (int j = 0; j < sampleSizesParallelOp[i]; j++) {
            cout << "checked" << endl;
            clearArray(matrixC_ParallelOp, (i + 1) * 100);
            averageParOp += multiplyMatrixParallelOp(sampleSizesParallelOp[i], (i + 1) * 100);
        }

        averageSer /= sampleSizesSerial[i];
        averagePar /= sampleSizesParallel[i];
        averageParOp /= sampleSizesParallelOp[i];

        cout << "Average time in sec to mulitply in Serial for " << (i + 1) * 100 << " inputs = " << averageSer << endl;
        cout << "Average time in sec to mulitply in Parallel for " << (i + 1) * 100 << " inputs = " << averagePar <<
        endl << endl;
        cout << "Average time in sec to mulitply in Optimized Parallel for " << (i + 1) * 100 << " inputs = " <<
        averageParOp << endl << endl;
    }


}

/*------------------------------------------------
        Optimized parallel multiplication
------------------------------------------------*/

void pMultiply(int n, double matrixA[1000][1000], double matrixB[1000][1000], double matrixC[1000][1000],
               double matrixT[1000][1000]) {
    int threadId = omp_get_thread_num();    // get thread id
    int len = n / 2;

    int ax_start = (threadId / 4) * len;    // starting x coordinate of sub matrix A
    int ay_start = ((threadId % 4) / 2) * len;    // starting y coordinate of sub matrix A
    int bx_start = (threadId % 2) * len;    // starting x coordinate of sub matrix B

    int i, j, k;

    int by_start = ax_start;    // starting y coordinate of sub matrix B
    int ay_end = ay_start + len;    // ending y coordinate of sub matrix A
    int bx_end = bx_start + len;    // ending x coordinate of sub matrix B

    /*multiply using matrix A and transpose of matrix B*/
        #pragma omp parallel for collapse(3)
        for (i = ay_start; i < ay_end; i++) {
            for (j = bx_start; j < bx_end; j++) {
                for (k = 0; k < len; k++) {
                    matrixC[i][j] += (matrixA[i][ax_start + k] * matrixB[j][by_start + k]);
                }
            }
        }
         #pragma omp parallel for collapse(3)//
        for (i = ay_start; i < ay_end; i++) {
            for (j = bx_start; j < bx_end; j++) {
                for (k = 0; k < len; k++) {
                    matrixT[i][j] += (matrixA[i][ax_start + k] * matrixB[j][by_start + k]);
                }
            }
        }
}

double multiplyMatrixParallelOp(int nSamples, int matrixSize) {
    int iteration = 0;
    int n = matrixSize;

    /*cout << "Number of samples: ";
    cin >> nSamples;*/

    double runningTime[nSamples];    // keep running times of each sample
    double average = 0;    // average value of running times

    int i, j, k;

    initializeMatrix(n);

    double start_time = omp_get_wtime();

    /*parallel for loop to get the transpose of matrix B*/
#pragma omp parallel for collapse(2)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matrixTransB[j][i] = matrixB[i][j];
        }
    }

    //Multiply matrices using parallel divide and conquer approach
    // Number of threads is 8

    pMultiply(n, matrixA, matrixTransB, matrixC_ParallelOp, matrixT);

    //Get resulting matrix for A*B using parallel for loop
#pragma omp parallel for collapse(2)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matrixC_ParallelOp[i][j] += matrixT[i][j];
        }
    }

    double end_time = omp_get_wtime();
    double timeTaken = end_time - start_time;

    //Clear allocated memory
    clearArray(matrixA, n);
    clearArray(matrixB, n);
    clearArray(matrixC_ParallelOp, n);
    clearArray(matrixT, n);
    clearArray(matrixTransB, n);
    return timeTaken;
}

/*function initialization*/



void initializeMatrix(int n) {
    //initialize both matrices with random values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrixA[i][j] = rand();
            matrixB[i][j] = rand();
        }
    }
}

double multiplyMatrixParallel(int n) {
    //matrix multiplication
    double begin = omp_get_wtime();    //start time
#pragma omp parallel for schedule(dynamic)		//this for loop will be evenly distributed among n number of threads threads and be calculated independent
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                matrixC_Parallel[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
    double end = omp_get_wtime();;        //end time

    return (double) (end - begin);        //duration for multiplication is calculated
}

double multiplyMatrixSerial(int n) {
    //matrix multiplication
    double begin = omp_get_wtime();    //start time
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                matrixC_Serial[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
    double end = omp_get_wtime();    //end time

    return (double) (end - begin);        //duration for multiplication is calculated
}

void clearArray(double matrix[1000][1000], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                matrix[i][j] = 0.0;
            }
        }
    }
}

