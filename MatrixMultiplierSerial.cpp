#include <iostream>
#include <cstdlib> 
#include <string>
#include <ctime>
#include <fstream>

using namespace std;

/*matrix declaration*/
double matrixA[1000][1000];
double matrixB[1000][1000];
double matrixC[1000][1000];
/*matrix declaration*/

void initializeMatrix(int n);
void writeMatrix(double matrix[1000][1000], char const* fileName ,int n);
void multiplyMatrix(int n);

int main(){
	int n = 0;
	cout << "Enter matrix dimension \n";
	cin >> n;

	initializeMatrix(n);
	cout << "Matrix A and Matrix B was initialized: \n";

	cout << "Matrix A Matrix B is being multiplied: \n";

	clock_t begin = clock();
	multiplyMatrix(n);
	clock_t end = clock();

	float duration = (float)(end-begin)/CLOCKS_PER_SEC;

	cout << "Time taken for multiplication : " << duration<<endl<<endl;

	cout << "Initialized Matrix A was written to file matrixA.txt: \n";
	writeMatrix(matrixA,"matrixA.txt",n);
	cout << "Initialized Matrix B was written to file matrixA.txt: \n";
	writeMatrix(matrixB,"matrixB.txt",n);
	cout << "Multiplied answer, Matrix C was written to file matrixC.txt: \n";
	writeMatrix(matrixC,"matrixC.txt",n);

}

void initializeMatrix(int n){
	
	for (int i=0; i<n; i++){		
		for (int j=0; j<n; j++){
			matrixA[i][j] = rand();
			matrixB[i][j] = rand();
		}
	}
}

void writeMatrix(double matrix[1000][1000], char const* fileName, int n){

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

void multiplyMatrix(int n){
	
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
				double sum = 0;
			for (int k = 0; k < n; k++){
				sum += matrixA[i][k] * matrixB[k][j];
			}		
			matrixC[i][j] = sum;
		}
		
	}
}