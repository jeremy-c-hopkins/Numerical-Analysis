#include <iostream>
#include <cmath>
using namespace std;

int main(){

    int n;
    int xEigen;
    int yEigen;
    int iter;

    cout << "Enter the node count (n^2): " << endl;
    cin >> n;
    cout << "Enter x eigen value: " << endl;
    cin >> xEigen;
    cout << "Enter y eigen value: " << endl;
    cin >> yEigen;
    cout << "Enter iteration count: " << endl;
    cin >> iter;

    double arr[n][n];
    double f[n];
    double h=M_PI/(n+1.0);
    int sqrtN=sqrt(n);

    // Setting up Laplace
    for(int i=0; i<n; i++){
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
        }
        f[i]=sin(xEigen*(i+1)*h)*sin(yEigen*(i+1)*h);
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(isnan(arr[i][j])){
                arr[i][j]=0;
            }
            cout << arr[i][j] << "\t";
        }
        cout << f[i] << "\n";
    }

    for(int i=0; i<n; i++){
        double sum = 0;
        for(int j=0; j<n; j++){
            
        }
    }
}
