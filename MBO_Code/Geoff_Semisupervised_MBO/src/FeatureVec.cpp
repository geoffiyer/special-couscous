// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#include "FeatureVec.h"
#include <cmath>
#include <iostream>
using namespace std;

FeatureVec::FeatureVec(){};
FeatureVec::~FeatureVec(){};

//input is a m by n by 3 matrix (3m by n), output give each frame a margin
//with f rows/columns. The output is a M by N by 3 matrix (3M by N).
//The (i,j,k)th element of input is input(k*m*n+i*N+j)
void FeatureVec::padarray(double* input, double* output, int f, int m, int n){
    int M = m + 2 * f;
    int N = n + 2 * f;
    
    for (int k=0; k<3; k++) {
        
        for (int j=f; j<N-f; j++) {
            for (int i=f; i<M-f; i++) {
                output[k*M*N+N*i+j] = input[k*m*n+n*(i-f)+(j-f)];
            }
        }
        
        for (int j=0; j<f; j++) {
            for (int i=f; i<M-f; i++) {
                output[k*M*N+N*i+j] = input[k*m*n+(i-f)*n+(f-j-1)];
                output[k*M*N+N*i+j+f+n] = input[k*m*n+(i-f)*n+(n-j-1)];
            }
        }
        
        for (int j=0; j<N; j++) {
            for (int i=0; i<f; i++) {
                output[k*M*N+i*N+j] = output[k*M*N+(2*f-i-1)*N+j];
                output[k*M*N+(i+f+m)*N+j] = output[k*M*N+(f+m-i-1)*N+j];
            }
        }
    }
}
// This function compute the gaussian kernel
// kernel is a 2f+1 by 2f+1 matrix
void FeatureVec::makeKernel (double* kernel, int f) {
    int m = 2*f+1;
    int n = 2*f+1;
    for (int d=1; d<f+1; d++) {
        double value = 1.0 / ((2.0*d+1.0)*(2.0*d+1.0));
        for (int i=-d; i<d+1; i++) {
            for (int j=-d; j<d+1; j++) {
                kernel[(f-i)*n+(f-j)] = kernel[(f-i)*n+(f-j)] + value;
            }
        }
    }

    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            kernel[i*n+j] /= f;
        }
    }
    
}
//This function builds the feature vector for each pixel in an RGB image
// input M(m+2f) by N(n+2f) by 3 matrix (3M by N), output a (m*n*3) by
//(2f+1)^2 feature matrix. Each row is a feature vector corrsponding to
//one pixel in the original m by n image
void FeatureVec::makeFeatureVec (double* Feature, double* paddedmatrix,
                                 double* kernel, int f, int m, int n) {
    //kernel is (2f+1)*(2f+1) by 1
    //paddedmatrix is (m+2f) by (n+2f) by 3 (3(m+2f) by (n+2f))
    //step1 feature vec based on neighbours
    // m*n by (2f+1)^2 by 3
    double* Feature_1 = new double[m*n*(2*f+1)*(2*f+1)*3];
    int M = m*n;
    int N = (2*f+1)*(2*f+1);
    for (int k=0; k<3; k++) {
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                int d = 0;
                for (int p = i; p < 2*f+1 +i; p++) {
                    for (int q = j; q < 2*f+1 +j; q++) {
                        Feature_1[k*M*N+(i*n+j)*N+d] =
                            paddedmatrix[k*(m+2*f)*(n+2*f)+p*(2*f+n)+q];
                        d += 1;
                    }
                }
            }
        }
    }
    // apply kernel
    for (int k=0; k<3; k++) {
        for (int i=0; i<M; i++) {
            for (int j=0; j<N; j++) {
                Feature_1[k*M*N+i*N+j] =
                    Feature_1[k*M*N+i*N+j] * sqrt(kernel[j]);
            }
        }
    }
    // reshape to (m*n) by ((2f+1)^2*3)
    for (int k=0; k<3; k++) {
        for (int i=0; i<M; i++) {
            for (int j=0; j<N; j++) {
                Feature[i*3*N+j+k*N] = Feature_1[k*M*N+i*N+j];
            }
        }
    }
    delete [] Feature_1;
}
