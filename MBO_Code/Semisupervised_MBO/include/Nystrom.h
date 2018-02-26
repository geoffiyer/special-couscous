// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#ifndef ____Nystrom6__
#define ____Nystrom6__



class Nystrom
{
public:
    Nystrom();
    ~Nystrom();
    void inverse(double* A, int N);
    void svd_solver(double * a,double * u, double * s, double * vt, int M);
    void Nystrom_algorithm(double * super_test, int m, double sigma, int N,
                           int d,double * V, double * D);

};

#endif /* defined(____Nystrom6__) */
