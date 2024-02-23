#include <iostream>
#include <cmath>
using namespace std;



int main(){
    int n;

    cout << "Enter matrix size n >= 3:" << endl;
    cin >> n;

    double arr[n][n];
    double lArr[n][n];
    double uArr[n][n];
    double solArr[n];
    double sol[n];
    double error = 0;

    // Hilbert matrix with sol x_4 = -2 and x_8 = 1
    for(int i=0; i<n; i++){
	    for(int j=0; j<n; j++){
	        arr[i][j]=1.0/(i+j+1.0);
	        lArr[i][j]=1.0/(i+j+1.0);
	    }
	    solArr[i]=1.0/(i+7+1.0)-2.0/(i+3+1.0);
	    sol[i]=0;
	    if(i==3){
	        sol[i]=-2.0;
	    }
	    if(i==7){
	        sol[i]=1.0;
	    }
    }

    // Printing Hilbert matrix
    cout << "Array" << "\n" <<endl;
    for(auto&j: lArr){
	    for(auto&i:j){
	        cout << i << "\t";
	    }
	    cout << "\n";
    }

    // Starting Cholesky decomposition
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += lArr[j][k] * lArr[j][k];
                }
                lArr[j][j] = sqrt(arr[j][j] - sum);
		        uArr[j][j] = lArr[j][j];
            }
            else {
                for (int k = 0; k < j; k++) {
                    sum += lArr[i][k] * lArr[j][k];
                }
                lArr[i][j] = (arr[i][j] - sum) / lArr[j][j];
		        uArr[j][i] = lArr[i][j];
            }
        }
	    for(int j=i+1; j<n; j++){
	        lArr[i][j]=0;
	        uArr[j][i]=0;
	    }
    }
    // Forward substitution
    for(int i=0; i<n; i++){
        solArr[i] = solArr[i]/lArr[i][i];
        for(int j=i+1; j<n; j++){
            lArr[j][i] = lArr[j][i] * solArr[i];
            solArr[j] -= lArr[j][i];
        }
    }
    // Backward substitution
    for(int i=n-1; i>-1; i--){
        solArr[i]=solArr[i]/uArr[i][i];
        for(int j=i-1; j>-1; j--){
            uArr[j][i] = uArr[j][i] * solArr[i];
            solArr[j] -= uArr[j][i];
        }
    }
    // Printing solution array
    for(int i=0; i<n; i++){
        cout << "x_" << i << ": " << solArr[i] << endl;
        error += pow(solArr[i] - sol[i], 2);
    }

    cout << "Cholesky Error: " << sqrt(error) << endl;
}


