// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#ifndef MBO_S_h
#define MBO_S_h
class MBO_S
{
public:
    MBO_S();
    ~MBO_S();
    void get_index_test(double *Index, double *f,int rows,int cols);
    void quickSort(double * arr, int left, int right);
    void projection(double * A, int rows, int cols);
    void initialization(double * u0, double * fidelity, int rows, int cols,
                        int rows_fidelity);
    void MBO_algorithm(double * index, double * V_test, double * D_test,
                       double * fidelity, double * u0, int N, int n,
                       int rows_fidelity, int M,int s, double dt, double c1,
                       double iterNum);
};

#endif /* MBO_S_h */
