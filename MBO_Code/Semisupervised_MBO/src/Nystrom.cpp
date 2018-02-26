// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.


#include <time.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Nystrom.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <cblas.h>
#include <omp.h>
using namespace std;

Nystrom::Nystrom(){};
Nystrom::~Nystrom(){};

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    
    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK,
                 int* lwork, int* INFO);
}

// A is an N by N matrix, calculate the inverse of A
void Nystrom::inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;
    
    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
    
    delete IPIV;
    delete WORK;
}

extern "C" void dgesdd_( char* jobz, int* m, int* n, double* a, int* lda,
                        double* s, double* u, int* ldu, double* vt, int* ldvt,
                        double* work, int* lwork, int* iwork, int* info );

//SVD a = vt*s*u
void Nystrom::svd_solver(double * a,double * u, double * s, double * vt, int M)
{
    int LDA=M;
    int LDU=M, LDVT=M, info, lwork;
    double wkopt;
    double* work;
    int* iwork= new int[8*M];
    int m = M, n = M, lda = LDA, ldu = LDU, ldvt = LDVT;
    lwork=-1;
    char jobz[]= "All";
    dgesdd_(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info );
    lwork = (int)wkopt;
    work = new double[lwork];
    
    dgesdd_(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info );
    
    if( info > 0 ) {
        cout<<"The algorithm computing SVD failed to converge.\n";
    }
    
    delete[] work;
    delete[] iwork;
}

extern "C"{
    double ddot_( const int *, const double *, const int *, const double *,
                 const int * );
}

void Nystrom::Nystrom_algorithm(double * super_test, int m, double sigma,
                                int N, int d,double * V, double * D)
{
    std::vector<int> myvector;
    for (int i=0; i<N; ++i) myvector.push_back(i);
    std::random_shuffle ( myvector.begin(), myvector.end() );
    
    // creating Wxx (m by m) and Wxy (m by N-m)
    
    double * Wxx = new double [m*m];
    int inca = 1, incb = 1;
    double num_1, num_2, num_3, dd;
    for(int i = 0; i<m; i++){
        int temp_1 = myvector[i]*d;
        num_2 = ddot_(&d, super_test+temp_1, &inca, super_test+temp_1, & incb);
        for(int j = 0; j<m; j++){
            int temp_2 = myvector[j]*d;
            num_1 = ddot_(&d, super_test+temp_1, &inca, super_test+temp_2, & incb);
            num_3 = ddot_(&d, super_test+temp_2, &inca, super_test+temp_2, & incb);
            dd = 1-num_1/sqrt(num_2*num_3);
            Wxx[i*m+j]=exp(-dd/sigma);
        }
    }
    // Wxy , m by N-m
    double * Wxy = new double [m*(N-m)];
    int j;
    
#ifdef OMP
    
#pragma omp parallel for default(shared) private(j,num_1,num_2,num_3,dd)  schedule(static)
    
    for(j = 0; j<N-m; j++){
        if(j==0){
            int N = omp_get_num_threads();
            cout<<"number of threads:"<<N<<endl;
        }
        
        int temp_2 = myvector[j+m]*d;
        num_3 = ddot_(&d, super_test+temp_2, &inca, super_test+temp_2, & incb);
        
        for(int i = 0; i<m; i++){
            int temp_1 = myvector[i]*d;
            
            num_1 = ddot_(&d, super_test+temp_1, &inca, super_test+temp_2, & incb);
            num_2 = ddot_(&d, super_test+temp_1, &inca, super_test+temp_1, & incb);
            
            dd = 1-num_1/sqrt(num_2*num_3);
            Wxy[i*(N-m)+j]=exp(-dd/sigma);
        }
    }
#else
    for(j = 0; j<N-m; j++){
        int temp_2 = myvector[j+m]*d;
        num_3 = ddot_(&d, super_test+temp_2, &inca, super_test+temp_2, & incb);
        
        for(int i = 0; i<m; i++){
            int temp_1 = myvector[i]*d;
            
            num_1 = ddot_(&d, super_test+temp_1, &inca, super_test+temp_2, & incb);
            num_2 = ddot_(&d, super_test+temp_1, &inca, super_test+temp_1, & incb);
            
            dd = 1-num_1/sqrt(num_2*num_3);
            Wxy[i*(N-m)+j]=exp(-dd/sigma);
        }
    }
#endif
    
    //dx = Wxx*1_m + Wxy * 1_{N-m}
    double * one_1 = new double [m];
    for(int i = 0; i<m; i++){ one_1[i] = 1; }
    double * one_2 = new double [N-m];
    for(int i = 0; i<N-m; i++) {one_2[i] = 1;}
    
    double * dx1 = new double [m];
    double * dx2 = new double [m];
    
    //dx1 = Wxx*one_1   m by 1
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    int incy = 1;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, 1, m, alpha, Wxx, m, one_1, 1, beta, dx1, 1);
    
    // dx2 = Wxy * 1_{N-m}    m by 1
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, 1, N-m, alpha, Wxy, N-m, one_2, 1, beta, dx2, 1);
    
    // dx2 = alpha*dx1+dx2
    
    cblas_daxpy(m,alpha,dx1,incx,dx2,incy);
    
    delete [] dx1;
    
    // keep the original copy of Wxx into Wxx0    
    double * Wxx0 = new double [m*m];
    
    int mm = m*m;
    int incOne = 1;
    cblas_dcopy(mm,Wxx,incOne,Wxx0,incOne);
    
    //Wxx0 = Wxx0^{-1}
    inverse(Wxx0,m);
    
    //dy = W_yx*W_xx^-1*W_xy*1_N-m
    //dy1 = Wxy*1_N-m, m by 1
    double * dy1 = new double [m];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, 1, N-m, alpha, Wxy, N-m, one_2, 1, beta, dy1, 1);
    delete [] one_2;
    //dy2 = Wxx0*dy1, m by 1
    double * dy2 = new double [m];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, 1, m, alpha, Wxx0, m, dy1, 1, beta, dy2, 1);
    delete [] Wxx0;
    delete [] dy1;
    // dy3 = Wxy^T*dy2, N-m by 1
    double * dy3 = new double [N-m];
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                N-m, 1, m, alpha, Wxy, N-m, dy2, 1, beta, dy3, 1);
    delete [] dy2;
    
    // dy4 = Wxy^T * 1_m
    double * dy4 = new double [N-m];
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                (N-m), 1, m, alpha, Wxy, (N-m), one_1, 1, beta, dy4, 1);
    
    delete [] one_1;
    
    // dy4 = alpha*dy3 + dy4, N-m by 1
    
    int n2 = N-m;
    cblas_daxpy(n2,alpha,dy3,incx,dy4,incy );
    
    delete [] dy3;
    
    // Normalizing Wxx and Wxy
    
    for(int i= 0; i<m; i++){
        dx2[i] = 1.0/sqrt(dx2[i]);
    }
    
    for(int i = 0; i<N-m; i++){
        dy4[i] = 1.0/sqrt(dy4[i]);
    }
    
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
            Wxx[i*m+j] = Wxx[i*m+j]*dx2[i]*dx2[j];
        }
    }
    
    for(int i=0; i<m; i++){
        for(int j=0; j<N-m; j++){
            Wxy[i*(N-m)+j] = Wxy[i*(N-m)+j]*dx2[i]*dy4[j];
        }
    }
    
    delete [] dx2;
    delete [] dy4;
    
    double * u = new double [m*m];
    double * s = new double [m];
    double * vt = new double [m*m];
    // copy Wxx into Wxx1 for solver, since Wxx1 would be changed by solver
    double * Wxx1 = new double [m*m];
    cblas_dcopy(mm,Wxx,incOne,Wxx1,incOne);
    
    svd_solver(Wxx1,u,s,vt,m);
    
    delete [] Wxx1;
    //u = Bx_t, vt = Bx, since u is stored column wise
    //Wxx = vt*s*u;
    //and Wxx1 is changed
    
    // s1 = s^-0.5, s2 = s^0.5
    double * s1 = new double [m];
    double * s2 = new double [m];
    
    for(int i = 0;i<m;i++){
        s1[i] = 1.0/sqrt(s[i]);
    }
    for(int i = 0;i<m;i++){
        s2[i] = sqrt(s[i]);
    }
    
    delete [] s;
    
    //S = vt*s*u
    // S1 = vt*diagnoal(s1)
    double * S1 = new double [m*m];
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
            S1[i*m+j] = vt[i*m+j]*s1[j];
        }
    }
    delete [] s1;
    
    // S = S1*u
    double * S = new double [m*m];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, m,m, alpha, S1, m, u, m, beta, S, m);
    
    // Q = Wxx +S(Wxy*Wxy^T)S
    // Q2 = Wxy*Wxy^T, m by m
    double * Q2 = new double [m*m];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                m, m, N-m, alpha, Wxy, N-m, Wxy, N-m, beta, Q2, m);
    
    // Q1 = S*Q2, m by m
    double * Q1 = new double [m*m];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, m, m, alpha, S, m, Q2, m, beta, Q1, m);
    
    delete [] Q2;
    // Q3 = Q1*S, m by m
    double * Q3 = new double [m*m];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, m, m, alpha, Q1, m, S, m, beta, Q3, m);
    delete [] Q1;
    delete [] S;
    
    cblas_daxpy(mm,alpha,Wxx,incx,Q3,incy );
    
    delete [] Wxx;
    
    double * R = new double [m*m];
    double * Rt = new double [m*m];
    double * Ksi = new double [m];
    // Q3 = Rt*Ksi*R
    svd_solver(Q3,R,Ksi,Rt,m);
    
    delete [] R;
    delete [] Q3;
    // copyint Ksi into Ksi1 so to calculate Ksi^-0.5 afterwords.
    double * Ksi1 = new double [m];
    
    // Q3 was changed by solver,
    
    double * Uupper = new double [m*m];
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
            Uupper[i*m+j] = vt[i*m+j]*s2[j];
        }
    }
    delete [] vt;
    delete [] s2;
    double * Ulower = new double [(N-m)*m];
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                N-m, m, m, alpha, Wxy, N-m, S1, m, beta, Ulower, m);
    delete [] Wxy;
    delete [] S1;
    
    // combine Uupper and Ulower
    
    double * U = new double [N*m];
    copy(Uupper,Uupper+m*m,U);
    copy(Ulower,Ulower+(N-m)*m,U+m*m);
    
    delete [] Uupper;
    delete [] Ulower;
    
    // making ksi a diagonal matrix
    for(int i = 0; i<m; i++){
        Ksi1[i] = 1.0/sqrt(Ksi[i]);
    }
    
    // V1 = Rt*diagonal(Ksi1)
    double * V1 = new double [m*m];
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
            V1[i*m+j] = Rt[i*m+j]*Ksi1[j];
        }
    }
    
    delete [] Ksi1;
    delete [] Rt;
    //V2 = u * V1
    double * V2 = new double [m*m];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, m, m, alpha, u, m, V1, m, beta, V2, m);
    delete [] V1;
    delete [] u;
    
    // V3 = U*V2,   N by m    U Nby m, V2 m by m
    double * V3 = new double [N*m];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, m, m, alpha, U, m, V2, m, beta, V3, m);
    
    delete [] V2;
    delete [] U;
    // unshuffle V
    // V is N by m from input
    for(int i = 0; i<N; i++){
        int temp = myvector[i];
        for(int j = 0; j<m; j++){
            V[temp*m+j] = V3[i*m+j];
        }
    }
    
    delete [] V3;
    
    // D = 1-Ksi
    // D is m by 1 from input
    for(int i = 0; i<m; i++){
        D[i] = 1-Ksi[i];
    }
    
    delete [] Ksi;
}






