#include <iostream>
#include <cstdlib> 
#include <string>
#include <ctime>
#include <fstream>
#include <omp.h>

using namespace std;

/*matrix declaration*/
int sampleSizes[10];
double matrixA[1000][1000];
double matrixB[1000][1000];
double matrixC[1000][1000];
double matrixD[1000][1000];
double matrixE[1000][1000];
/*matrix declaration*/

/*function initialization*/
void initializeMatrix(int n);
void writeMatrix(double matrix[1000][1000], char const* fileName ,int n);
double multiplyMatrixSerial(int n);
double multiplyMatrixParallel(int n);
/*function initialization*/

int main(){
	int n = 0;
	cout << "Enter matrix dimension \n";
	cin >> n;	// less than 1000, greater than 0

	initializeMatrix(n);	//initialize matrices with random values
	cout << "Matrix A and Matrix B was initialized: \n";

	cout << "Matrix A Matrix B is being multiplied serially: \n";

	
	double durationSerial = multiplyMatrixSerial(n);		//multiply matrices serially, which returns time taken
	
	cout << "Time taken for multiplication : " << durationSerial << endl << endl;

	cout << "Matrix A Matrix B is being multiplied serially: \n";

	
	double durationParallel = multiplyMatrixParallel(n);		//multiply matrices parallely, which returns time taken
	
	cout << "Time taken for multiplication : " << durationParallel << endl << endl;

	//write initialized and multiplied matrices to files
	cout << "Initialized Matrix A was written to file matrixA.txt: \n";
	writeMatrix(matrixA,"matrixA.txt",n);
	cout << "Initialized Matrix B was written to file matrixA.txt: \n";
	writeMatrix(matrixB,"matrixB.txt",n);
	cout << "Multiplied answer, Matrix C was written to file matrixC.txt: \n";
	writeMatrix(matrixC,"matrixC.txt",n);
	cout << "Multiplied answer, Matrix D was written to file matrixC.txt: \n";
	writeMatrix(matrixD,"matrixD.txt",n);

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

void writeMatrix(double matrix[1000][1000], char const* fileName, int n){
	//write matrix elements to files
	ofstream myfile (fileName); 
	if(myfile.is_open()){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){			
				myfile << matrix[i][j] << "\t" ;
			}	
			myfile << endl;
		}
		myfile.close();
	}
	else{
		cout <<"Unable to open "<< fileName <<endl;
	}
	
}

double multiplyMatrixParallel(int n){
	//matrix multiplication
	double begin = omp_get_wtime();	//start time
	#pragma omp parallel for collapse(3) schedule(dynamic)		//this for loop will be evenly distributed among n number of threads threads and be calculated independent
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
