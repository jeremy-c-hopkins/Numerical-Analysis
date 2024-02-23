#include <iostream>
#include <cmath>
using namespace std;

/*
 * Solves Hilbert matrix with and without partial pivoting.
 *
 * @version February 5, 2024
 * @author Jeremy Hopkins
 */

int main() {
    int n;
    int pivot;
    double max;
    double temp;
    double error = 0;

    cout << "Matrix n:" << endl; cin>>n;
    cout << "0 - no pivoting, 1 - standard pivoting, 2 - scaled partial pivoting:" << endl; cin>>pivot;

    double arr[n][n];
    double solArr[n];
    double scale[n];
    double solutionSpace[n];

    // Test matrix Hilbert 13 and 7th minus 2 times 3rd. - @author Robert Palais
    for(int i=0; i < n;i++){
        for(int j=0;j<n;j++){
            arr[i][j]=1.0/(i+j+1.0);
        }
        solArr[i]=1.0/(i+7+1.0)-2.0/(i+3+1.0);
        scale[i]=1.0/(i+7+1.0)-2.0/(i+3+1.0);
	    solutionSpace[i] = 0;
	    if(i == 3){
	        solutionSpace[i] = -2;
	    }
	    if(i == 7){
	        solutionSpace[i] = 1;
	    }
    }

    // Displaying Hilbert matrix. - @author Robert Palais
    cout << "Matrix and right side:" << "\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout << arr[i][j] <<"\t";
        }
        cout <<"\t \t" << solArr[i] << "\n"; } cout << "\n";

    // Gauss reduction with and without scaled partial pivoting.
    for (int i = 0; i < n; i++){

        // Checking for pivot conditions.
        if (pivot == 1){
            int index = -1;
            max = arr[i][i];
            for (int k = i + 1; k < n; k++) {
                // Checking if values below diagonal are larger than diagonal.
                if (arr[k][i] > max) {
                    max = abs(arr[k][i]);
                    index = k;
                }
            }
            // Checking if no values below diagonal were larger than diagonal
            if (index != -1){
                for (int k = i; k < n; k++){
                    swap(arr[i][k], arr[index][k]);
                }
                swap(solArr[i], solArr[index]);
                cout << "swap " << i+1 << " and " << index+1 << endl;
            }
        }

        if (pivot == 2){
            max = arr[i][i]/scale[i];
            int indexVal = 0;
            for (int j = i + 1; j < n; j++){
                temp = arr[j][i]/scale[j];
                if (temp > max){
                    max = temp;
                    indexVal = j;
                    for (int k = 0; k < n; k++){
                        swap(arr[j][k], arr[i][k]);
                        swap(scale[j], scale[i]);
                        swap(solArr[j], solArr[i]);
                    }
                }
            }
	    cout << "Swap " << i+1 << " and " << indexVal+1 << endl;
        }

        // Doing reduction for column i.
        for (int j = i + 1; j < n; j++){
            if (arr[j][i] == 0){
                continue;
            }
            double mij = arr[j][i] / arr[i][i];
            solArr[j] -= mij * solArr[i];
            for (int k = i; k < n; k++){
                arr[j][k] -= mij * arr[i][k];
            }
        }
    }

    // Back substitution
    for(int i = n-1; i > -1; i--){
        solArr[i] = solArr[i] / arr[i][i];
        for(int j = i - 1; j > -1; j--){
            arr[j][i] = arr[j][i] * solArr[i];
            solArr[j] -= arr[j][i];
        }
    }


    cout << "Solution array" << endl;

    // Printing solution array.
    for (int x = 0; x < n; x++){
        cout << "x_" << x << ": " << solArr[x] <<endl;
        error += pow(solutionSpace[x] - solArr[x], 2);
    }

    cout << "Error: " << sqrt(error) << endl;
}
