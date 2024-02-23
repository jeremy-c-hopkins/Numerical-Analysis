#include <iostream>
#include <cmath>
#include <iomanip> //for setting precision
using namespace std;

// g++ jacobi_gauss-seidel.cpp  -lm -o  jacobi
// ./jacobi

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
    int sqrtN=sqrt(n);
    double h=M_PI/(sqrtN+1.0);
    

    // Setting up Laplace
    for(int i=0; i<n; i++){
        for (int l=0;l<n;l++) {
 		    arr[i][l]=0.0;
 		} 
        arr[i][i] = -4.0;

        
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
            if(j+sqrt(n)==i){
                arr[i][j]=1.0;
            }
            if(abs(arr[i][j])<0.1){
                arr[i][j]=0;
            }
            if (i<sqrtN and j<sqrtN){
                int index = i*sqrtN +j;
                f[index]=sin(xEigen*(i+1)*h)*sin(yEigen*(j+1)*h);
            }
        }
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(isnan(arr[i][j])){
                arr[i][j]=0;
            }
            cout << arr[i][j] << "  ";
        }
        cout << "\t" << setprecision(5)<<f[i] << "\n";
        cout << "\t" << "\n";
    }

    for(int i=0; i<n; i++){
        double sum = 0;
        for(int j=0; j<n; j++){
            
        }
    }
}
