#include <iostream>
#include <cstdlib>
#include <string>
#include <ctime>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

/*hold calculated sample sizes CI 95% MoE 5% for mini samples*/
int sampleSizesParallelOp[10];

//to calculate sample size, average and sd is needed
const int sampleSizeEst = 10;

double parallelTimeStampOp[sampleSizeEst];


/*matrix declaration*/
double matrixA[1000][1000];
double matrixB[1000][1000];
double matrixC_ParallelOp[1000][1000];
double matrixD[1000][1000];
double matrixE[1000][1000];
double matrixT[1000][1000];
double matrixTransB[1000][1000];

/*function initialization*/
void initializeMatrix(int n);
void clearArray(double matrix[1000][1000], int n);
double multiplyMatrixParallelOp(int nSamples, int matrixSize);
void pMultiply(int n, double matrixA[1000][1000], double matrixB[1000][1000], double matrixC[1000][1000],
               double matrixT[1000][1000]);

int main() {

    initializeMatrix(1000);    //initialize matrices with random values
    cout << "Matrix A and Matrix B was initialized: \n";

    for (int i = 0; i < 10; i++) {
        double averageParOp = 0.0;
        double sdParOp = 0.0;

        for (int j = 0; j < sampleSizeEst; j++) {
            clearArray(matrixC_ParallelOp, (i + 1) * 100);
            parallelTimeStampOp[j] = multiplyMatrixParallelOp(sampleSizeEst, (i + 1) * 100);
            averageParOp += parallelTimeStampOp[j];
        }

        for (int j = 0; j < sampleSizeEst; j++) {
            sdParOp += pow((parallelTimeStampOp[j] - averageParOp), 2);
        }

        sdParOp = sqrt(sdParOp);

        sampleSizesParallelOp[i] = (int) pow(((100.0 * 1.96 * (double) sdParOp) / (5 * averageParOp)), 2);

        sampleSizesParallelOp[i] = max(sampleSizesParallelOp[i], sampleSizeEst);


        cout << "Optimized Parallel" << "SD " << sdParOp << " AVG " << averageParOp << endl;

        cout << "Sample size for mulitply in Optimized Parallel " << (i + 1) * 100 << " = " <<
        sampleSizesParallelOp[i] << endl << endl;
    }

    cout << "Calculating sample size is finished" << endl;

    for (int i = 0; i < 10; i++) {
        clearArray(matrixC_ParallelOp, (i + 1) * 100);

        double averageSer = 0.0;

        double averagePar = 0.0;

        double averageParOp = 0.0;

        for (int j = 0; j < sampleSizesParallelOp[i]; j++) {
            cout << "checked" << endl;
            clearArray(matrixC_ParallelOp, (i + 1) * 100);
            averageParOp += multiplyMatrixParallelOp(sampleSizesParallelOp[i], (i + 1) * 100);
        }

        averageParOp /= sampleSizesParallelOp[i];

        cout << "Average time in sec to mulitply in Parallel for " << (i + 1) * 100 << " inputs = " << averagePar <<
        endl << endl;
        cout << "Average time in sec to mulitply in Optimized Parallel for " << (i + 1) * 100 << " inputs = " <<
        averageParOp << endl << endl;
    }


}

void initializeMatrix(int n) {
    //initialize both matrices with random values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrixA[i][j] = rand();
            matrixB[i][j] = rand();
        }
    }
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


