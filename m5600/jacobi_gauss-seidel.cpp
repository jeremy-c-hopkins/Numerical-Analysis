#include <iostream>
#include <cmath>
#include <iomanip> //for setting precision
using namespace std;

// g++ jacobi_gauss-seidel.cpp  -lm -o  jacobi
// ./jacobi

// This will print out Laplace matrix and eigenfunction on its right side.
// On next line, it will print datas with appropriate column on top.


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
    double newGuessArr[n];
    double oldGuessArr[n];
    double gsArr[n];
    int sqrtN=sqrt(n);
    double h=M_PI/(sqrtN+1.0);
    double error = 100;
    

    // Setting up Laplace
    for(int i=0; i<n; i++){
        newGuessArr[i]=0; //initialize
        oldGuessArr[i]=0; //initialize
        gsArr[i]=0; //initialize

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

    // Printing Laplace matrix and solution
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

    
    // Random
    for (int k=0; k<iter;k++){ //iteration
        //Jacobi
        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<n; j++){
                if(j != i){
                    sum += arr[i][j] * newGuessArr[j];
                }
            }
            newGuessArr[i] = (f[i] - sum)/arr[i][i];
            error = oldGuessArr[i] - newGuessArr[i];
            oldGuessArr[i] = newGuessArr[i];
        }
        //Gauss Seidel
        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<n; j++){
                if(j != i){
                    sum += arr[i][j] * gsArr[j];
                }
            }
            gsArr[i] = (f[i] - sum)/arr[i][i];
        }
    }
    
    //print
    cout<<"jacobi      |   jacobi over eigenvector y    |   Gauss-Seidel    | Gauss Seidel over eiven vector y"<<"\n";
    for(int i=0; i<n; i++){
        cout<<newGuessArr[i]<<"\t"<<"\t"<<newGuessArr[i]/f[i]<<"\t"<<"\t"<<gsArr[i]<<"\t"<<"\t"<<gsArr[i]/f[i]<<"\n";
    }
    cout << "\n" << error << endl;

}
