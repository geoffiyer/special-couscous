// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.


#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "MBO_S.h"
#include <cblas.h>
using namespace std;

// #define PRINT_MBO

MBO_S::MBO_S(){};
MBO_S::~MBO_S(){};
//f is a rows by cols matrix, index is a rows by 1 matrix,
//we find the largest element in each row of f and put its index in Index.
void MBO_S::get_index_test(double *Index, double *f,int rows,int cols) {
    for(int i = 0; i<rows; i++) {
        Index[i] = 0;
        double max = f[i*cols];
        for(int j = 0; j<cols; j++) {
            if(max<f[i*cols+j]) {
                max = f[i*cols+j];
                Index[i] = j;
            }
        }
    }
}

void MBO_S::MBO_algorithm(double * index, double * V, double * D, double *
                    fidelity, double * u0, int N, int n, int rows_fidelity,
                    int M,int s, double dt, double C, double iterNum) {

    /* Geoff's notes to self:
       index: class label for each pixel. Updated throughout loop
       V: The eigenvectors. I don't know the format
       D: The eigenvalues. I don't know the format
       fidelity: the fidelity
       u0: the starting class assigment. Mostly random except for the
           semisupervised part.
       N: number of points
       n: Number of classes
       rows_fidelity: number of fidelity points
       M: Number of eigenvectors
       s: Number of diffusions between each threshold
       dt: stepsize
       C: same as mu from main.cpp. The fidelity term constant
       iterNum: max number of iterations
    */

#ifdef PRINT_MBO
    ofstream mout;
    mout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/MBOclass.txt");
    if(mout.fail()) {
        cout << "MBO output failed to open\n";
    }
#endif

    // initialization
    // For each row of u0, we randomly select one component to be 1 and others
    // to be 0. For the rows corresponding to fidelity points, we assign the
    // component corresponding to the known labels to be 1 and others 0.
    for(int i = 0; i<N*n; i++){
        u0[i] = 0;
    }    for(int i = 0; i<N; i++)
    {
        int randint = rand()%n;
        int x = i*n+randint;
        u0[x] = 1;
    }
    for(int i = 0; i<rows_fidelity; i++){
        for(int j = 0; j<n; j++){
            int temp = fidelity[2*i]*n+j;
            u0[temp]=0;
        }
        int temp2 = fidelity[2*i]*n+fidelity[2*i+1];
        u0[temp2] = 1;
    }

    // ////////////////////////////////////////////////
    // // Print initial values just to do a check:
    // ofstream fout;
    // fout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt");
    // if(fout.fail()) {
    //     cout << "MBO initialization output failed to open\n";
    // }
    // get_index_test(index, u0,N,n);
    // for(int i=0; i<N; ++i) {
    //     fout << index[i] << "\n";
    // }
    // // End print operation
    // ///////////////////////////////////////////////
    
    // some constants
    double alpha = 1.0;
    double beta = 0.0;
    double gamma = -1.0;
    int incOne = 1;
    double theta = -dt/s;
    int len_u = N*n;
    int len_a = M*n;
    
    // allocating memory
    double* a = new double[M*n];
    double* Denom = new double[M];
    double* d = new double[M*n];
    double* u = new double[N*n];
    double* y0 = new double [N*n];
    double* X = new double[N*n];
    double* u_old = new double[N*n];
    double* u_diff = new double[N*n];
    
    // in implicit formula, a_k^{n+1} = 1/(1+dt*lambda_k+dt*C)*a_k^n-dt*d_k^n,
    // Denom[i] = 1/(1+dt*lambda_k+dt*C)
    for(int i =0; i<M; i++){
        Denom[i] = 1.0/(1.0+D[i]*dt/s+C*dt/s);
    }
    
    // initialize the coefficients C*lambda(x)*(u-\hat{u}) = d*V
    for(int i= 0; i<M*n; i++){
        d[i] = 0;
    }
    // u = u0
    cblas_dcopy(len_u,u0,incOne,u,incOne);
    
    // u_old = u0
    cblas_dcopy(len_u,u0,incOne,u_old,incOne);
    
#ifdef PRINT_MBO
    // mout << "C: " << C << "\n";
    // mout << "s: " << s << "\n";
    // mout << "len_u: " << len_u << "\n";
    // mout << "len_a: " << len_a << "\n";
    mout << "The Denom:\n";
    for(int i=0; i<M; ++i)
        mout << Denom[i] << "\n";
    // mout << "\nHere's the starting u:\n";
    // for(int i=0; i<N; ++i) {
    //     for(int j=0; j<n; ++j)
    //         mout << u[j+i*n] << " ";
    //     mout << "\n";
    // }
    mout << endl;
#endif

    ///////////iteration  ////////////////////
    for(int iter=0; iter<iterNum;iter++){
        // a is the coefficient of u projected on the eigenvectors V,
        // u = V*a, so a = V^T*u
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    M, n, N, alpha, V, M, u, n, beta, a, n);
        // step 1: diffuse
        // we do step1 (diffuse) s times before step 2 (thresholding)
        for(int j = 0; j<s;j++){
            // a = Denom*(a-(dt/s)*d)
            cblas_daxpy(len_a,theta,d,incOne,a,incOne );
            for(int i = 0; i<M; i++){
                for(int j = 0; j<n;j++){
                    a[i*n+j]=a[i*n+j]*Denom[i];
                }
            }
            
            // u = V*a;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        N, n, M, alpha, V, M, a, n, beta, u, n);
            // copy u into y0, y0 is used for u = u-u0, u is for the
            // thresholding
            cblas_dcopy(len_u,u,incOne,y0,incOne);
            
            //y0 = u-u0
            cblas_daxpy(len_u,gamma,u0,incOne,y0,incOne );

            // In the end, X = u-lambda(x)*(u-u0). lambda(x) = 1 if this pixel has
            // known label, and 0 otherwise
            // Part 1: Start with X = 0
            for(int i = 0; i<N*n; i++){
                X[i]=0;
            }
            // Part 2: X = - lambda(x)(u-u0})
            for(int i = 0; i<rows_fidelity; i++){
                for(int j = 0; j<n; j++){
                    int temp = fidelity[i*2]*n+j;
                    X[temp] = -y0[temp];
                }
            }
            // Part 3: X_final = u + X
            cblas_daxpy(len_u,alpha,u,incOne,X,incOne);
            
            // d is the coefficient of C*(u-lambda(x)*(u-u0)) projected onto V,
            // C*X = V*d, so d = C*V^T*X
            // Part 1: d = V^T*X
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                M, n, N, alpha, V, M, X, n, beta, d, n);
            // Part 2: d = d*C;
            for(int i = 0; i<M*n; i++){
                d[i] = d[i]*C;
            }
        }

#ifdef PRINT_MBO
        mout << "u before threshold:\n";
        for(int i=0; i<N; ++i) {
            for(int j=0; j<n; ++j)
                mout << u[j+i*n] << " ";
            mout << "\n";
        }
        mout << endl;
#endif

        //step 2: thresholding part
        get_index_test(index, u,N,n);
        
        for(int i = 0; i<N*n; i++) {
            u[i]=0;
        }
        
        for(int i = 0; i<N; i++) {
            int x = i*n+index[i];
            u[x]= 1;
        }

        //// criteria to stop
        // norm1 = ||u||_2^2
        double norm1 = cblas_dnrm2(len_u,u,incOne);
        
        // u_diff = u
        cblas_dcopy(len_u,u,incOne,u_diff,incOne);
        
// #ifdef PRINT_MBO
//         mout << "u after threshold:\n";
//         for(int i=0; i<N; ++i) {
//             for(int j=0; j<n; ++j)
//                 mout << u[j+i*n] << " ";
//             mout << "\n";
//         }
//         mout << "norm is " << norm1 << "\n";
//         mout << endl;
// #endif

        // u_diff = u_diff-u_old
        cblas_daxpy(len_u,gamma,u_old,incOne,u_diff,incOne );
        
        // norm2 = ||u-u_old||_2^2
        double norm2 = cblas_dnrm2(len_u,u_diff,incOne);
        
        if(norm2/norm1<0.0000001){
            cout << "Num iter: " << iter << "\n";
            break;
        }
        //u_old = u;
        cblas_dcopy(len_u,u,incOne,u_old,incOne);

    }

    delete [] a;
    delete [] Denom;
    delete [] d;
    delete [] u;
    delete [] y0;
    delete [] X;
    delete [] u_old;
    delete [] u_diff;
    
}
