// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#ifndef FeatureVec_h
#define FeatureVec_h

class FeatureVec {
public:
    FeatureVec();
    ~FeatureVec();
    void padarray (double* input, double* output, int f, int m, int n);
    void makeKernel (double*kernel, int f);
    void makeFeatureVec (double* Feature, double* paddedmatrix,
                            double* kernel, int f, int m, int n);
};
#endif /* FeatureVec_h */
