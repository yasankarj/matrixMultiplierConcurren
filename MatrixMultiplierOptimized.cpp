#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

using namespace std;

/*Get random double values in range min and max exclusive*/
double GetRandom(double min, double max) {
    double f = (double) rand() / RAND_MAX;
    return min + f * (max - min);
}

/*Parallel mutipication single thread function*/
void pMultiply(int n, double **matrixA, double **matrixB, double **matrixC, double **matrixT) {
    int threadId = omp_get_thread_num();    // get thread id
    int len = n / 2;

    /*Calculate resulting subsection of the matrices A,B and C
    for given thread using thread id*/

    int ax_start = (threadId / 4) * len;    // starting x coordinate of sub matrix A


    int ay_start = ((threadId % 4) / 2) * len;    // starting y coordinate of sub matrix A


    int bx_start = (threadId % 2) * len;    // starting x coordinate of sub matrix B

    int i, j, k;

    int by_start = ax_start;    // starting y coordinate of sub matrix B
    int ay_end = ay_start + len;    // ending y coordinate of sub matrix A
    int bx_end = bx_start + len;    // ending x coordinate of sub matrix B

    if(threadId/4==0){// first 4 threads calculate matrix C values
        for (i = ay_start; i < ay_end; i++) {
            for (j = bx_start; j < bx_end; j++) {
                for (k = 0; k < len; k++) {
                    /*calculate matrix multipication using matrix A
                    and transpose of matrix B*/
                    matrixC[i][j] += (matrixA[i][ax_start + k] * matrixB[j][by_start + k]);
                }
            }
        }
    }
    else{
        for (i = ay_start; i < ay_end; i++) {
            for (j = bx_start; j < bx_end; j++) {
                for (k = 0; k < len; k++) {
                    /*calculate matrix multipication using matrix A
                    and transpose of matrix B*/
                    matrixT[i][j] += (matrixA[i][ax_start + k] * matrixB[j][by_start + k]);
                }
            }
        }
    }

}

void pMultiplyForT(int n, double **matrixA, double **matrixB, double **matrixC, double **matrixT) {

    int threadId = omp_get_thread_num();    // get thread id
    //cout<<"thread"<<threadId<<endl;
    int len = n / 2;

    /*Calculate resulting subsection of the matrices A,B and C
    for given thread using thread id*/
    //int bin_ax1 = (threadId/4);
    int ax_start = (threadId / 4) * len;    // starting x coordinate of sub matrix A

    //int bin_ay1 = (threadId%4)/ 2;
    int ay_start = ((threadId % 4) / 2) * len;    // starting y coordinate of sub matrix A

    //int bin_bx1 = (threadId%2);
    int bx_start = (threadId % 2) * len;    // starting x coordinate of sub matrix B

    int i, j, k;

    int by_start = ax_start;    // starting y coordinate of sub matrix B
    int ay_end = ay_start + len;    // ending y coordinate of sub matrix A
    int bx_end = bx_start + len;    // ending x coordinate of sub matrix B

    for (i = ay_start; i < ay_end; i++) {
        for (j = bx_start; j < bx_end; j++) {
            for (k = 0; k < len; k++) {
                /*calculate matrix multipication using matrix A
                and transpose of matrix B*/
                matrixT[i][j] += (matrixA[i][ax_start + k] * matrixB[j][by_start + k]);
            }
        }
    }

}

int run(int nSamples, int matrixSize, bool isTrial) {

    //int nSamples = numberOfSamples;
    int iteration = 0;
    int n = matrixSize;
    bool checked = false;
    int requiredNumberOfSamples = 30;

    /*cout << "Number of samples: ";
    cin >> nSamples;*/

    double runningTime[nSamples];    // array to keep running times of each sample
    double avarage = 0;    // average value of running times

    int i, j, k;

    /*Declare double 2D arrays for matrices*/
    double **matrixA = new double *[n];
    double **matrixB = new double *[n];
    double **matrixC = new double *[n];
    double **matrixT = new double *[n];
    double **matrixTransB = new double *[n];

    for (i = 0; i < n; i++) {
        matrixA[i] = new double[n];
        matrixB[i] = new double[n];
        matrixC[i] = new double[n];
        matrixT[i] = new double[n];
        matrixTransB[i] = new double[n];
    }

    /*Iterarte through each sample*/
    for (iteration = 0; iteration < nSamples; iteration++) {

        srand(time(NULL));    // seed with time to generate different random values

        /*Initialize A and B matrices with random values*/
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                matrixA[i][j] = GetRandom(0, 10);
                matrixB[i][j] = GetRandom(0, 10);
            }
        }

        //struct timeval startTime, endTime, resultTime;

      //  gettimeofday(&startTime, NULL);    // get starting time of calculations
        double start_time = omp_get_wtime();

        /*use parallel for loop to get the transpose of matrix B*/
#pragma omp parallel for collapse(2) private (i,j) schedule(dynamic)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                matrixTransB[j][i] = matrixB[i][j];
            }
        }

        /*Compute matrix multipication using parallel divide and conquer method
        with 8 parallel threads */
#pragma omp parallel num_threads(8)
        pMultiply(n, matrixA, matrixTransB, matrixC, matrixT);
/*#pragma omp parallel num_threads(4)
        pMultiplyForT(n, matrixA, matrixTransB, matrixC, matrixT);*/

        /*calculate final matirix answer using parallel for loop*/
#pragma omp parallel for collapse(2) private (i,j) schedule(dynamic)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                matrixC[i][j] += matrixT[i][j];
            }
        }

        double end_time = omp_get_wtime();
       // gettimeofday(&endTime, NULL);    // get end time of calculations

        /*resultTime.tv_sec = endTime.tv_sec - startTime.tv_sec;
        resultTime.tv_usec = endTime.tv_usec - startTime.tv_usec;*/
        double timeTaken = end_time - start_time;
        //cout<<timeTaken<<endl;

        /*double totalTime = (double) resultTime.tv_sec * 1000 +
                           (double) resultTime.tv_usec / 1000; */   // calculate resulting time in miliseconds

        // runningTime[iteration] = totalTime;
        runningTime[iteration] = timeTaken;
        //avarage += totalTime;
        avarage+=timeTaken;
    }

    double sd = 0;

    //cout<<"average :"<<avarage<<endl;
    //cout<<"samples :"<<nSamples<<endl;
    avarage = avarage / nSamples;    // calculate the mean

    /*Calculate standard deviation*/
    for (i = 0; i < nSamples; i++) {
        sd += pow((runningTime[i] - avarage), 2);
    }

    sd = sd / nSamples;
    sd = sqrt(sd);

    if(isTrial){
        cout << "----------------------"<<endl;
        cout << "For sample size of 30,"<<endl;
    }
    else{
        cout << "----------------------"<<endl;
        cout << "For required sample size of "<<nSamples <<endl;
    }
    cout << "	Mean: " << avarage << endl;
    cout << "	SD: " << sd << endl;
    cout << "	Average Running Time : "<<avarage<<endl;

/*Calculate required samples to hold 5% accuracy with 95% confidence intervel*/
    if (isTrial) {
        requiredNumberOfSamples = (int) pow(((100.0 * 1.96 * sd) / (5 * avarage)), 2);
        cout << "Required Samples: " << requiredNumberOfSamples << endl;
        cout << "----------------------"<<endl;
    }


/*Clear allocated memory*/
    for (i = 0; i < n; i++) {
        delete[] matrixA[i];
        delete[] matrixB[i];
        delete[] matrixC[i];
        delete[] matrixT[i];
        delete[] matrixTransB[i];
    }
    delete[] matrixA;
    delete[] matrixB;
    delete[] matrixC;
    delete[] matrixT;
    delete[] matrixTransB;

    return requiredNumberOfSamples;
}

int main() {

    int nSamples = 30;
    int iteration = 0;
    int n = 0;
    bool checked = false;

    cout << "Matrix size: ";
    cin >> n;

    nSamples = run(30,n,true);
    run(nSamples,n,false);
}


