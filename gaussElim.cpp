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
    int count = 0;
    double max;
    double temp;
    double error = 0;

    cout << "Matrix n:" << endl; cin>>n;
    cout << "0 - no pivoting, 1 - standard pivoting, 2 - scaled partial pivoting:" << endl; cin>>pivot;

    double arr[n][n];
    double arrInverse[n][n];
    double l[n][n];
    double lInverse[n][n];
    double tempSolArr[n];
    double solArr[n];
    double scale[n];
    double solutionSpace[n];

    // Test matrix Hilbert 13 and 7th minus 2 times 3rd. - @author Robert Palais
    for(int i=0; i < n;i++){
        for(int j=0;j<n;j++){
            arr[i][j]=1.0/(i+j+1.0);
            if(i!=j){
                l[i][j] = 0;
                lInverse[i][j] = 0;
                arrInverse[i][j] = 0;
            }
        }
        l[i][i] = 1.0;
        lInverse[i][i] = 1.0;
        arrInverse[i][i] = 1.0;
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
        cout <<"\t \t" << solArr[i] << "\n"; 
    } cout << "\n";

    // Gauss reduction with and without scaled partial pivoting.
    if(pivot==1||pivot==2){
        for (int i = 0; i < n; i++){

            // // Checking for pivot conditions.
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

            // Scaled partial pivoting
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

            for(int i=0; i<n; i++){
                for(int i=0; i<n; i++){
                    // Doing reduction for column i.
                    for (int j = i + 1; j < n; j++){
                        if (arr[j][i] == 0){
                            continue;
                        }
                        double mij = arr[j][i] / arr[i][i];
                        solArr[j] -= mij * solArr[i];
                        for (int k = i; k < n; k++){
                            arr[j][k] -= mij * arr[i][k];
                            count++;
                        }
                    }   
                }
            }
        }

        // Back substitution
        for(int i = n-1; i > -1; i--){
            solArr[i] = solArr[i] / arr[i][i];
            for(int j = i - 1; j > -1; j--){
                arr[j][i] = arr[j][i] * solArr[i];
                solArr[j] -= arr[j][i];
                count++;
            }
        }

    } else{
        for(int i=0; i<n; i++){
            // Doing reduction for column i.
            for (int j = i + 1; j < n; j++){
                double mij = arr[j][i] / arr[i][i];
                l[j][i] = mij;
                for (int k = i; k < n; k++){
                    arr[j][k] -= mij * arr[i][k];
                    count++;
                }
            }   
        }

        cout << "\nL MATRIX\n";
        for(auto&i:l){
            for(auto&j:i){
                cout << j << " ";
            }
            cout << "\n";
        }
        cout << "\n";

        cout << "U MATRIX\n";
        for(auto&i:arr){
            for(auto&j:i){
                cout << j << " ";
            }
            cout << "\n";
        }
    
        // Inverting the Lower Matrix
        for(int i = 1; i<n; i++){
            for(int j = 0; j < i; j++){
                double mij = l[i][j];
                for(int k = 0; k < i; k++){
                    lInverse[i][k] -= mij * lInverse[j][k];
                    l[i][k] -= mij * l[j][k];
                    count++;
                }
            }
        }

        // Scaling the values on the upper matrix diagonal to 1
        for (int i = n-1; i>-1; i--){
            arrInverse[i][i] = 1/arr[i][i];
            double scale = 1/arr[i][i];
            for(int l=i; l<n; l++){
                arr[i][l] = scale * arr[i][l];
                count++;
            }
        }

        // Inverting the upper matrix
        for (int i = 0; i<n-1; i++){
             for (int j = i+1; j<n; j++){
                double mij = arr[i][j];
                for (int k = j; k<n; k++){
                    arr[i][k] -= mij * arr[j][k];
                    arrInverse[i][k] -= mij * arrInverse[j][k];
                    count++;
                }
            }
        }

        cout << "\n\nL INVERSE\n";
        for(auto&i:lInverse){
            for(auto&j:i){
                cout << j << "  ";
            }
            cout << "\n";
        }
        cout << "\nU INVERSE\n";
        for(auto&i:solArr){
            cout << i << "\n";
        }
        cout << "\n";


        

        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<n; j++){
                sum += lInverse[i][j] * solArr[j];
                count++;
            }
            tempSolArr[i] = sum;
        }

        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<n; j++){
                sum += arrInverse[i][j] * tempSolArr[j];
                count++;
            }
            solArr[i] = sum;
        }
        
    }


    
    //=======================================================================
    //=======================================================================
    cout << "Solution" << endl;

    // Printing solution array.
    for (int x = 0; x < n; x++){
        cout << "x_" << x << ": " << solArr[x] <<endl;
        error += pow(solutionSpace[x] - solArr[x], 2);
    }

    cout << "Error (for n>=7): " << sqrt(error) << endl;
    cout << "Number of Operations: " << count;
}
