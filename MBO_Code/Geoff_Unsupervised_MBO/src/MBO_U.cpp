// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#include "MBO_U.h"
#include <iomanip>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
extern "C" {
#include <cblas.h>
}

#define PRINT_MBO

using namespace std;

MBO_U::MBO_U(){};
MBO_U::~MBO_U(){};

// x and c has same number of cols. result is a matrix with dimension x_row
// by c_row. the i,j element of result is ||x_i-c_j||^2, which is the
// Euclidean distance of the ith row of x and the jth row of c.
void dist2(double* x, double* c, double* result, int &x_row, int &x_col,
           int &c_row, int &c_col){
    int incx = 1;
    int len = x_row*c_row;
    double alpha = 1.0;
    double theta = -2.0;
    double beta = 0.0;
    // ||x||^2
    double * tempx = new double[len];
    
    for(int i=0; i<x_row; i++){
        double norm = cblas_dnrm2(x_col,x+i*x_col,incx);
        for(int j=0; j<c_row; j++){
            tempx[i*c_row+j] = norm;
        }
    }
    // x*y^T
    double * tempxy = new double[x_row*c_row];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, x_row, c_row,
                x_col, alpha, x, x_col, c, x_col, beta, tempxy, c_row);
    
    //  ||c||^2
    for(int i=0; i<c_row; i++){
        double norm2 = cblas_dnrm2(c_col,c+i*c_col,incx);
        for(int j=0; j<x_row; j++){
            result[j*c_row+i] = norm2;
        }
    }
    cblas_daxpy(len,alpha,tempx,incx,result,incx);
    cblas_daxpy(len,theta,tempxy,incx,result,incx);
}

//f is a rows by cols matrix, index is a rows by 1 matrix, we find
//the largest element in each row of f and put its index in Index.
void MBO_U::get_index(double *Index, double *f,int rows,int cols){
    for(int i = 0; i<rows; i++){
        Index[i] = 0;
        double max = f[i*cols];
        for(int j = 0; j<cols; j++){
            if(max<f[i*cols+j]){
                max = f[i*cols+j];
                Index[i] = j;
            }
        }
    }
}
// calculate the purity score between two vectors
double MBO_U::puritymeas(double * index, double * preindex, int len,
                         int num_class){
    double score = 0;
    double avgscore = 0;
    for(int i = 0; i< num_class; i++){
        int pos = 0;
        std::vector<double> ind;
        for(int j = 0; j<len; j++){
            if(index[j] == i){
                ind.push_back(pos);
            }
            pos++;
        }
        int count = 0;
        int vec_length = ind.size();
        for(int k = 0; k<vec_length; k++){
            int temp = ind.at(k);
            if(preindex[temp]==i)
                count++;
        }
        score = score + count;
    }
    avgscore = score/len;
    return avgscore;
}

void MBO_U::Mumford_Shah(double * final_index, double * super_test,
                         double * V_test, double * D_test, int N, int n,
                         int M, int d, double dt, double lambda){
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    
    // initialization
    double * ff_test = new double [N*n];
    double *index_test = new double[N];
    double * preindex_test = new double [N];
    for(int i = 0; i<N*n; i++){
        ff_test[i] = 0;
    }
    for(int i = 0; i<N; i++){
        int randint = rand()%n;
        index_test[i] = randint;
        int x = i*n+randint;
        ff_test[x] = 1;
    }

#ifdef PRINT_MBO
    ofstream mout;
    mout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/MBOclass.txt");
    if(mout.fail()) {
        cout << "MBO output failed to open\n";
    }
    #endif
    
    cblas_dcopy(N,index_test,incx,preindex_test,incx);
    //a = V_test^T * ff_test  M by n matrix
    double * a_test =new double[M*n];
    
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                M, n, N, alpha, V_test, M, ff_test, n, beta, a_test, n);

    double * Denom_test = new double[M];
    for(int i = 0; i<M ; i++){
        Denom_test[i] = 1.0 - dt*D_test[i];
    }
    
    delete [] D_test;
    
    int flag = 1;
    int iter = 1;
    
    double * f1 = new double[N*n];
    double * c_test = new double[n*d];
    double * ff_sum_test = new double[n];
    double * cr = new double[d];
    double * f2 = new double[N*n];

    /* Geoff's notes to self:
       ff_test is current classification (zeros and ones)
               starts out random then get's edited in each MBO loop
               in the paper this is called u^n
       a_test  is the same as the a from the paper. It's like ff but in evec coords
       c_test  is the average value of the original input data in each class
               (ex: row 1 of c_test is  average(input data that is currently in class 1)
               Why do we use this??
       f1      is the intermediate classification (probability. decimal values)
               They call this u^{n+ 1/2}
       f2      is made up of differences between the input data value and
               the c_test (average per class) value
               THIS IS WHAT LAMBDA IS USED FOR
    */
    
    while(flag){
        for(int i = 0; i<M; i++){
            for(int j = 0; j<n; j++){
                a_test[i*n+j] = a_test[i*n+j]*Denom_test[i];
            }
        }

        //f1 = V_test*a_test
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N, n, M, alpha, V_test, M, a_test, n, beta, f1, n);
        
        // updating c
        // c_test = ff_test^T * super_test
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,n, d, N,
                    alpha, ff_test, n, super_test, d, beta, c_test, d);
        
        // number of nodes in each class <1,ff_test>
        for(int i=0; i<n; i++){
            ff_sum_test[i] = 0;
        }
        for(int i=0; i<N; i++){
            int idx = index_test[i];
            ff_sum_test[idx] += 1;
        }

#ifdef PRINT_MBO
        // Printout:
        mout << setprecision(2) << fixed;
//    cout << "Printing MBO iter " << iter << "\n";
        for(int i=0; i<n; ++i) {
//            for(int j=0; j<n; j++) {
                mout << ff_sum_test[i] << " ";
//            }
            mout << "\n";
        }
        mout << "\n\n\n";
#endif
        
        // c_test is n by d, ff_sum_test is 1 by n/ n by 1,
        // Divide the ith row of c_test by the ith element of ff_sum_test
        for(int i = 0; i<n; i++){
            for(int j = 0; j<d; j++){
                c_test[i*d+j] = c_test[i*d+j]/ff_sum_test[i];
            }
        }
        
        dist2(super_test,c_test,f2,N,d,n,d);
        
        //f1 = -dt*lambda*f2+f1;
        int f_len = N*n;
        double gamma = -dt*lambda;
        // f1 = gamma*f2+f1
        cblas_daxpy(f_len,gamma,f2,incx,f1,incx);

        get_index(index_test,f1,N,n);
        
        for(int i = 0; i<N*n; i++){
            ff_test[i] = 0;
        }
        for(int i = 0; i<N; i++){
            int x = i*n+index_test[i];
            ff_test[x]= 1;
        }
        
        // a = V^T*ff
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,M, n, N,
                    alpha, V_test, M, ff_test, n, beta, a_test, n);
        
        //stopping creteria
        double purity_score_test = puritymeas(index_test,
                                              preindex_test, N, n);
        iter = iter + 1;
        
        if(purity_score_test>0.999){
            flag = 0;
        }
        
        cblas_dcopy(N,index_test,incx,preindex_test,incx);
    }

    cout << "Number of MBO iterations: " << iter-1 << "\n";
    
    final_index[0]= N;
    final_index[1]= 1;
    for(int i = 0; i<N; i++){
        final_index[i+2] = index_test[i];
    }
    
    //deallocate
    delete [] index_test;
    delete [] preindex_test;
    delete [] ff_test;
    delete [] a_test;
    delete [] Denom_test;
    delete [] f1;
    delete [] c_test;
    delete [] ff_sum_test;
    delete [] cr;
    delete [] f2;
}
