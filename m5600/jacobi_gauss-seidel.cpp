#include <iostream>
#include <cmath>
#include <iomanip> //for setting precision
using namespace std;

// g++ jacobi_gauss-seidel.cpp  -lm -o  jacobi
// ./jacobi

// This will print out Laplace matrix and eigenfunction on its right side.
// On next line, it will print 


int main(){
    int n; //matrix size
    int xEigen;
    int yEigen;
    int iter;

    cout << "Enter the node count (which equals to n^2): " << endl;
    cin >> n;
    cout << "Enter x eigen value: " << endl;
    cin >> xEigen;
    cout << "Enter y eigen value: " << endl;
    cin >> yEigen;
    cout << "Enter iteration count: " << endl;
    cin >> iter;

    double arr[n][n];
    double f[n];
    double guessArr[n][n];
    int sqrtN=sqrt(n);
    double h=M_PI/(sqrtN+1.0);
    double error = 100;
    

    // Setting up Laplace
    for(int i=0; i<n; i++){
        for (int l=0;l<n;l++) {
 		    arr[i][l]=0.0; //initialize
 		} 
        arr[i][i] = -4.0; //initialize

        if(i-1>-1 && i%sqrtN!=0){
            arr[i][i-1]=1.0;
        }
        if(i+1<n && ((i+1)%sqrtN)!=0){
            arr[i][i+1]=1.0;
        }

        for(int j=0; j<n; j++){
            if(j-sqrt(n)==i){
                arr[i][j]=1.0;
            }
            else if(j+sqrt(n)==i){
                arr[i][j]=1.0;
            }
            else if(abs(arr[i][j])<0.1){
                arr[i][j]=0;
            }
            // setting eigenfunction
            if (i<sqrtN and j<sqrtN){
                int index = i*sqrtN +j;
                f[index]=sin(xEigen*(i+1)*h)*sin(yEigen*(j+1)*h); 
            }
        }
    }

    


    //functions to print onto console/terminal
    cout<<"Laplace matrix  |   eigenfunction"<<"\n";
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(isnan(arr[i][j])){
                arr[i][j]=0;
            }
            cout << arr[i][j] << "  ";
        }
        cout << "\t" << setprecision(6)<<f[i] << "\n";
        cout << "\t" << "\n";
    }

    while(error > 1){
        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<i; j++){
                sum += arr[i][j] * guessArr[j];
            }
            for(int j=i+1; j<n; j++){
                sum += arr[i][j] * guessArr[j];
            }
            guessArr = (f[j] - sum)/(-4);
        }
    }
}
