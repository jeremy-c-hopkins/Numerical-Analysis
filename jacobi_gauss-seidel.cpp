#include <iostream>
#include <cmath>
#include <iomanip> //for setting precision
using namespace std;

// g++ jacobi_gauss-seidel.cpp  -lm -o  jacobi
// ./jacobi

/**
 * This will compute the Laplacian matrix and solve for Ax=y, with
 * y being sin(xh)sin(yh).
 * 
*/


int main(){
    int n; //matrix size
    int xEigen;
    int yEigen;
    int iter;

    cout << "\n\n";
    cout << "Enter node count of 4 to see single iteration of Ax=y with y being eigen vector of the Laplacian on a 2x2 grid." << "\n" <<endl;

    cout << "Enter the node count (n^2; i.e 4, 9, 16...): " << endl;

    cin >> n;
    int sqrtN=sqrt(n);

    cout << "Enter x eigen value from 1 to " << sqrtN << ": " << endl;
    cin >> xEigen;
    cout << "Enter y eigen value from 1 to " << sqrtN << ": " << endl;
    cin >> yEigen;
    cout << "Enter iteration count: " << endl;
    cin >> iter;

    
    // Setting up Ax and Y
    double A[n][n];
    double jacobiNewX[n];
    double jacobiOldX[n];
    double gaussSeidelX[n];
    double gaussSeidelX1[n];
    double gaussSeidelX2[n];
    double gaussSeidelX3[n];
    double gaussSeidelX4[n];
    double y[n];
    
    // Defining h value and setting desired error
    double h=M_PI/(sqrtN+1.0);

    double error = 100;


    /*
        Setting up Laplace with n^2 nodes. There will be 
        -4 across the diagonal and 1 at adjacent nodes for all n>=2.
    */ 
    for(int i=0; i<n; i++){
        jacobiNewX[i]=0; 
        jacobiOldX[i]=0; 
        gaussSeidelX[i]=0; 
        gaussSeidelX1[i]=0;
        gaussSeidelX2[i]=0;
        gaussSeidelX3[i]=0;
        gaussSeidelX4[i]=0;

        for (int l=0;l<n;l++) {
 		    A[i][l]=0.0; 
 		} 
        A[i][i] = -4.0; 

        if(i-1>-1 && i%sqrtN!=0){
            A[i][i-1]=1.0;
        }
        if(i+1<n && ((i+1)%sqrtN)!=0){
            A[i][i+1]=1.0;
        }

        for(int j=0; j<n; j++){
            if(j-sqrt(n)==i){
                A[i][j]=1.0;
            }
            else if(j+sqrt(n)==i){
                A[i][j]=1.0;
            }

            // Setting f(x)=y vector
            if (i<sqrtN and j<sqrtN){
                int index = i*sqrtN +j;
                y[index]=sin(xEigen*(i+1)*h)*sin(yEigen*(j+1)*h); 
            }
        }
    }

    // Printing Laplace matrix and solution
    cout << "\n" <<"Laplace Matrix   x   |   y" <<"\n";
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << A[i][j] << "  ";
        }
        cout << "  \tx_" << i+1 << " \t" << setprecision(2)<< y[i] << "\n";
        cout << "\t" << "\n";
    }

    
    for (int k=0; k<iter;k++){ // Iteration
        //Jacobi
        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<n; j++){
                if(j != i){
                    sum += A[i][j] * jacobiNewX[j];
                }
            }
            jacobiNewX[i] = (y[i] - sum)/A[i][i];
            error = min(error, jacobiOldX[i] - jacobiNewX[i]);
            jacobiOldX[i] = jacobiNewX[i];
        }

        //Gauss Seidel
        for(int i=0; i<n; i++){
            double sum = 0;
            for(int j=0; j<n; j++){
                if(j != i){
                    sum += A[i][j] * gaussSeidelX[j];
                }
            }
            gaussSeidelX[i] = (y[i] - sum)/A[i][i];
        }
    }


    // Print
    cout << "For " << iter << " iterations, the last iteration error is: " << error << "\n\n";
    cout<<"Jacobi      |   Jacobi over eigenvector y    |   Gauss-Seidel    | Gauss Seidel over eigen vector y"<<"\n";
    for(int i=0; i<n; i++){
        cout<< setprecision(3) <<jacobiNewX[i]<<"\t\t\t"<< setprecision(3) <<
        jacobiNewX[i]/y[i]<<" \t\t\t "<< setprecision(3) <<gaussSeidelX[i]<<"\t\t\t\t"<<
        setprecision(3) << gaussSeidelX[i]/y[i]<<"\n";
    }

    cout << "\n\n";

    /* 
        Showing one step iteration with eigen vector 
        of the Laplacian as solution.
    */ 
    if (n == 4){
        double eigenVectorY1[n] = {1, 1, 1, 1};
        double eigenVectorY2[n] = {1, -1, -1, 1};
        double eigenVectorY3[n] = {0, -1, 1, 0};
        double eigenVectorY4[n] = {-1, 0, 0, 1};

        double xEigenSolution1[n] = {-0.5, -0.5, -0.5, -0.5};
        double xEigenSolution2[n] = {-0.167, 0.167, 0.167, -0.167};
        double xEigenSolution3[n] = {0, 0.25, -0.25, 0};
        double xEigenSolution4[n] = {0.25, 0, 0, -0.25};


        for(int k=0; k<iter; k++){
            //Gauss Seidel
            for(int i=0; i<n; i++){
                double sum1 = 0;
                double sum2 = 0;
                double sum3 = 0;
                double sum4 = 0;

                for(int j=0; j<n; j++){
                    if(j != i){
                        sum1 += A[i][j] * gaussSeidelX1[j];
                        sum2 += A[i][j] * gaussSeidelX2[j];
                        sum3 += A[i][j] * gaussSeidelX3[j];
                        sum4 += A[i][j] * gaussSeidelX4[j];
                    }
                }
                gaussSeidelX1[i] = (eigenVectorY1[i] - sum1)/A[i][i];
                gaussSeidelX2[i] = (eigenVectorY2[i] - sum2)/A[i][i];
                gaussSeidelX3[i] = (eigenVectorY3[i] - sum3)/A[i][i];
                gaussSeidelX4[i] = (eigenVectorY4[i] - sum4)/A[i][i];
            }
        }

        cout << "\tSingle iteration for Ax = Y = eigen vector {1, 1, 1, 1}" << "\n";
        cout << "Single Iteration   |   Solution" << endl;
        for(int i=0; i<n; i++){
            cout << setprecision(3) << gaussSeidelX1[i] << "\t\t\t" << setprecision(3) << xEigenSolution1[i] << "\n";
        }
        cout << "\n\n";
        cout << "\tSingle iteration for Ax = Y = eigen vector {1, -1, -1, 1}" << "\n";
        cout << "Single Iteration   |   Solution" << endl;
        for(int i=0; i<n; i++){
            cout << setprecision(3) << gaussSeidelX2[i] << "\t\t\t" << setprecision(3) << xEigenSolution2[i] << "\n";
        }
        cout << "\n\n";
        cout << "\tSingle iteration for Ax = Y = eigen vector {0, -1, 1, 0}" << "\n";
        cout << "Single Iteration   |   Solution" << endl;
        for(int i=0; i<n; i++){
            cout << setprecision(3) << gaussSeidelX3[i] << "\t\t\t" << setprecision(3) << xEigenSolution3[i] << "\n";
        }
        cout << "\n\n";
        cout << "\tSingle iteration for Ax = Y = eigen vector {-1, 0, 0, 1}" << "\n";
        cout << "Single Iteration   |   Solution" << endl;
        for(int i=0; i<n; i++){
            cout << setprecision(3) << gaussSeidelX4[i] << "\t\t\t" << setprecision(3) << xEigenSolution4[i] << "\n";
        }
        cout << "\n\n";
        cout << "Checking Orthogonality\n";
        double innerproduct1 = 0;
        double innerproduct2 = 0;
        double innerproduct3 = 0;
        double innerproduct4 = 0;
        double innerproduct5 = 0;
        double innerproduct6 = 0;
        for(int i=0; i<n; i++){
            innerproduct1 += eigenVectorY1[i] * eigenVectorY2[i];
            innerproduct1 += eigenVectorY1[i] * eigenVectorY3[i];
            innerproduct1 += eigenVectorY1[i] * eigenVectorY4[i];
            innerproduct1 += eigenVectorY2[i] * eigenVectorY3[i];
            innerproduct1 += eigenVectorY2[i] * eigenVectorY4[i];
            innerproduct1 += eigenVectorY3[i] * eigenVectorY4[i];
        }
        cout << "Inner product of Eigenvector 1 and 2: " << innerproduct1 << "\n";
        cout << "Inner product of Eigenvector 1 and 3: " << innerproduct2 << "\n";
        cout << "Inner product of Eigenvector 1 and 4: " << innerproduct3 << "\n";
        cout << "Inner product of Eigenvector 2 and 3: " << innerproduct4 << "\n";
        cout << "Inner product of Eigenvector 2 and 4: " << innerproduct5 << "\n";
        cout << "Inner product of Eigenvector 3 and 5: " << innerproduct6 << "\n";
    }
}
