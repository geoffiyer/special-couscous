// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#ifndef ____Mumford_Shah6__
#define ____Mumford_Shah6__

class MBO_U
{
public:
    MBO_U();
    ~MBO_U();
    void get_index(double *Index, double *f,int rows,int cols);
    void Mumford_Shah(double * final_index, double * super_test,
                      double * V_test, double * D_test, int N, int n,
                      int M, int d, double dt, double lambda);
    double puritymeas(double * index, double * preindex, int len, int num_class);

};
#endif /* defined(____Mumford_Shah6__) */
