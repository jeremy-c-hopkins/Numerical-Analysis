#include <iostream>
#include <math.h>
using namespace std;
const int N=256;;

int main(){
double a[N][N],y[N],x[N],xgs[N],xoldj[N],xnewj[N],sum,tol=1e-15;
double pi=4.0*atan(1.0),h,hsq,p,q,f;
int i,j,k,l,n,nsq,count,iter,maxiter;
int component[N][4],neighbors[N];

cout << "Enter number of interior grid points per side of square, n:" << endl; cin>>n;
nsq=n*n; // matrix size
h=pi/(n+1.0); hsq = h*h; // grid spacing and square for denominator of approximation
cout << "Enter x eigenmode, from 1 to "<< n << ": " << endl; cin>>p;
cout << "Enter y eigenmode, from 1 to "<< n << ": " << endl; cin>>q;
cout << "Enter number of iterations:"  << endl; cin>>maxiter;

k=0; // counter for vector index
for(i=0;i<n;i++){ 
	for(j=0;j<n;j++){
	// for i,j+1 i,j-1  i+1,j i-1,j check if in interior and if so count and find index
	    count = 0;
		if (j+1<n){ 
            component[k][count]=k+1; 
            count++;
        }
		if (i+1<n){ 
            component[k][count]=k+n; 
            count++;
        }
		if (j-1>=0){ 
            component[k][count]=k-1; 
            count++;
        }
		if (i-1>=0){ 
            component[k][count]=k-n; 
            count++;
        }
		neighbors[k]=count;
		// construct right side: eigenvector sin(px)sin(qy);
		y[k]=sin(p*(i+1)*h)*sin(q*(j+1)*h);
		k++; 
        } 
    }
for (k=0;k<nsq;k++){
    xgs[k]=xoldj[k]=0.0;
}

cout << "\n";
cout << "F(x) values: " << "\n" << endl;
for(auto&i: y){
    cout << i << endl;
}
cout << "End" << endl;
cout << "\n" << endl;
// initialize Jacobi and Gauss-Seidel iterations
// This is not needed for algorithm!! Just for seeing the matrix!
for (k=0;k<nsq;k++) for (l=0;l<nsq;l++) a[k][l]=0.0; // initialize all zeros
for (k=0;k<nsq;k++) a[k][k]=-4.0; // initialize diagonal -4s
for (k=0;k<nsq;k++) for (l=0;l<neighbors[k];l++) a[k][component[k][l]]=1.0;
cout << "System for Laplacian on nxn square interior grid points " << "\n";
cout << "A: \n \n";
for(i=0;i<nsq;i++){ for(j=0;j<nsq;j++){
cout << a[i][j] <<"\t"; } cout << endl;} cout << "\n";

cout << "Right side eigenfunction:" << endl;
for (k=0;k<nsq;k++){
cout << y[k] << "\n";
}
// Iterations
for (iter=0;iter<maxiter;iter++){ // usually a while loop with convergence condition
	// Jacobi - don't update on the fly
	for (k=0;k<nsq;k++) {
		sum=0.0;
	    for (l=0;l<neighbors[k];l++)
            sum+=xoldj[component[k][l]];
		    xnewj[k]=-0.25*(hsq*y[k]-sum); 
        }
	for (k=0;k<nsq;k++){
        xoldj[k]=xnewj[k];


	// Gauss-Seidel - update on the fly
	for (k=0;k<nsq;k++) {
		sum=0.0;
		for (l=0;l<neighbors[k];l++){
            sum+=xgs[component[k][l]];
        }
		xgs[k]=-0.25*(hsq*y[k]-sum); 
    }
}

cout<< "Jacobi last iterate" << endl;

for (k=0;k<nsq;k++){
    cout<< xnewj[k]<< endl;
}

cout << endl;

cout<< "Gauss-Seidel last iterate" << endl;
for (k=0;k<nsq;k++){
cout<< xgs[k]<< endl;
}

cout<< "Jacobi Solution over eigenvector y" << endl;

for (k=0;k<nsq;k++){
cout<< xnewj[k]/y[k]<< endl;
}

cout << endl;

cout<< "Gauss-Seidel Solution over eigenvector y" << endl;

for (k=0;k<nsq;k++){
cout<< xgs[k]/y[k] << endl; 
}

return 0;
}
}
