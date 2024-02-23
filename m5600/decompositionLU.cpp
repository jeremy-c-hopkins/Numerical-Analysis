//
// Created by Jeremy Hopkins on 2/7/24.
//
#include "iostream"
using namespace std;


int main(){
    int n;
    int pivot;

    cout << "Enter matrix size n >=3:" << endl; cin>>n;
    cout << "Enter 0 for no pivoting, Enter 1 for pivoting:" << endl; cin>>pivot;

    double upperArr[n][n];
    double lowerArr[n][n];
    double solArr[n];

    // Test matrix Hilbert 13 and 7th minus 2 times 3rd. - @author Robert Palais
    for(int i=0; i < n;i++){
        lowerArr[i][i] = 1;
        for(int j=0;j<n;j++){
            upperArr[i][j]= 1.0 / (i + j + 1.0);
            if(j!=i){
                lowerArr[i][j] = 0;
            }
        }
        solArr[i]=1.0/(i+7+1.0)-2.0/(i+3+1.0);
    }

    cout << "Matrix and right side:" << "\n";
    for(int i=0;i<n;i++) {
        for (int j = 0; j < n; j++) {
            cout << upperArr[i][j] << "\t";
        }
        cout << "\n";
//        cout <<"\t \t" << solArr[i] << "\n";
    }
    cout << "\n";


    // LU Decomposition

    for (int i = 0; i < n; i++){

        // Doing reduction for column i.
        for (int j = i + 1; j < n; j++){
            if (upperArr[j][i] == 0){
                continue;
            }
            double mij = upperArr[j][i] / upperArr[i][i];
            solArr[j] -= mij * solArr[i];
            lowerArr[j][i] = mij;
            for (int k = i; k < n; k++){
                upperArr[j][k] -= mij * upperArr[i][k];
            }
        }
    }

    // Displaying Hilbert matrix. - @author Robert Palais
    cout << "Lower Matrix:" << "\n";
    for(int i=0;i<n;i++) {
        for (int j = 0; j < n; j++) {
            cout << lowerArr[i][j] << "\t";
        }
        cout << "\n";
//        cout <<"\t \t" << solArr[i] << "\n";
    }
    cout << "\n";


    cout << "Upper Matrix:" << "\n";
    for(int i=0;i<n;i++) {
        for (int j = 0; j < n; j++) {
            cout << upperArr[i][j] << "\t";
        }
        cout << "\n";
//        cout <<"\t \t" << solArr[i] << "\n";
    }
    cout << "\n";

    return 0;
}
