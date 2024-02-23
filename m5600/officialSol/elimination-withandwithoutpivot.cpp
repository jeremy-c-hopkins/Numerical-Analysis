#include <iostream>
#include <math.h>
using namespace std;
const int N = 40;

int main(){
int i,j,k, n, imax, tmpi,pivot; // loop indices and matrix size, output loops don't disturb
double dkk, mik, sum, tmp, max;  // multipliers reciprocals accumulator
int in[N];
double a[N][N],y[N],x[N],scale[N]; // matrix, data, solution entries, multipliers

        cout << "Enter matrix size n >=3:" << endl; cin>>n;
        cout << "Enter 0 for no pivoting, Enter 1 for pivoting:" << endl; cin>>pivot;
        // Test matrix Hilbert 13 and 7th minus 2 times 3rd
                for(i=0;i<n;i++){
                        for(j=0;j<n;j++){
                                a[i][j]=1.0/(i+j+1.0);
                                }
                        y[i]=1.0/(i+7+1.0)-2.0/(i+3+1.0);
                        }
    // display
        cout << "Matrix and right side:" << "\n";
        for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                        cout << a[i][j] <<"\t";
                }
                cout <<"\t \t" << y[i] << "\n"; } cout << "\n";

// Elimination with Implicit Scaled Partial Pivoting Algorithm Begins
        for(i=0;i<n;i++){ // determine scales
                max=0.0;
                for(j=0;j<n;j++){
                        tmp=fabs(a[i][j]);
                        if(tmp>max)max=tmp;
                        }
                // tmp=fabs(y[i]);if(tmp>max)max=tmp; // include RHS?
                scale[i]=max;
                in[i]=i;
        }
        for(k=0;k<n-1;k++){ // Scaled partial pivoting for column k
                // Implicit Scaled Partial Pivoting
                max=fabs(a[in[k]][k])/scale[in[k]];imax=k;
                for(i=k+1;i<n;i++){ // Pivoting
                        tmp=fabs(a[in[i]][k])/scale[in[i]];
                        if(tmp>max){max=tmp;imax=i; }
                }
                if (pivot !=0){
                if(imax!=k){tmpi=in[k];in[k]=in[imax];in[imax]=tmpi; cout<<"Swap "<<k<<","<<imax<<endl;
                }
        } // End pivot step for column k
                // Elimination for column k
                for(i=k+1;i<n;i++){
                        mik=a[in[i]][k]/a[in[k]][k];
                        for(j=k+1;j<n;j++) a[in[i]][j]-=mik*a[in[k]][j];
                        y[in[i]]-=mik*y[in[k]];a[in[i]][k]=mik;
                        }
                } // Elimination complete

                // Backsubstitution
                for(i=n-1;i>=0;i--){
                        sum=y[in[i]];
                        for(j=i+1;j<n;j++) sum-=a[in[i]][j]*x[j];
                        x[i]=sum/a[in[i]][i];
                        } // Solution complete
        cout << "Solution of Hilbert System:" << "\n";
        for(i=0;i<n;i++){ cout << x[i] <<"\n"; } cout << "\n";
return 0;
}
