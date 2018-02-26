// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include "MBO_S.h"
#include "Nystrom.h"
#include <string>
#include <vector>
#include "FeatureVec.h"
extern "C" {
#include "iio.h"
}
using namespace std;

// this function is for writting the result
// the first two elements of A are the dimentions
std::ostream& operator<<(std::ostream& out, double *A){
    int rows;
    int cols;
    rows = A[0];
    cols = A[1];
    
    for (int i=0;i<rows;i++)			{
        for (int j=0;j<cols;j++)
            out<<A[i*cols+j+2]<<"\t";
        out<<"\r\n";
    }
    out<<"\r\n\n";
    
    return out;
}

int main(int argc, char * argv[]){
    if (argc<6) {
        cout << "Not enough input parameters"<<endl;
        printf("Usage: %s n dt mu sigma input_image_string input_fidelity_string\n", argv[0]);
    } else {
        int n = atoi(argv[1]);
        double dt = atof(argv[2]);
        double mu = atof(argv[3])*100.0;
        double sigma = atof(argv[4]);
        string input_image_string(argv[5]);
        string input_fidelity_string;
    
        if (argc==7) {
            input_fidelity_string = argv[6];
        }
        double iterNum= 500;
        int s=3;
        int m = 100;
        int windowSize = 2;
        /////////////Input image///////////////
        int width;
        int height;
        int d;
        double * input_image;
        input_image = iio_read_image_double_vec(input_image_string.c_str(),
                                                &width,&height,&d);
        // add a small number to all the pixel values in case there are a lot of 0s
        for (int i=0; i<width*height*d; i++) {
            input_image[i] += 0.1;
        }
        double * feature_vec;
        
        if (d<4) {
        // for RGB image, rescale some parameters, so that it is
        // more consistent with the hyperspectral images
	        sigma = sigma*100.0;
            mu = mu*10.0;
        // for RGB image, need to read frame by frame and build the feature vector
            double * input_image_split = new double [width*height*d];
            
            for (int i=0; i<width*height; i++) {
                for (int l=0; l<d; l++) {
                    input_image_split[width*height*l+i] = input_image[d*i+l];
                }
            }
            
            FeatureVec Feature_obj;
            
            double * kernel = new double[(2*windowSize+1)*(2*windowSize+1)];
            for (int i=0; i<(2*windowSize+1)*(2*windowSize+1); i++) {
                kernel[i] = 0;
            }
            Feature_obj.makeKernel(kernel,windowSize);
            
            int N1p = height + 2 * windowSize;
            int N2p = width + 2 * windowSize;
            double* padimage = new double[N1p*N2p*3];
            for (int i=0; i<N1p*N2p*3; i++) {
                padimage[i] = 0;
            }
            Feature_obj.padarray(input_image_split,padimage,
                                 windowSize,height,width);
            
            int size_feature = width*height*(2*windowSize+1)*(2*windowSize+1)*3;
            feature_vec = new double[size_feature];
            
            Feature_obj.makeFeatureVec(feature_vec,padimage,kernel,
                                       windowSize,height,width);
            
            d = (2*windowSize+1)*(2*windowSize+1)*3;
            delete [] padimage;
        } else {
            cout << "Hyperspectral Image"<<endl;
            feature_vec = input_image;
        }
        
        int N = width*height;
        // input fidelity
        int width_f;
        int height_f;
        int d_f;
        vector<double> fidelity_vec;

        if (argc < 7) {
        } else {
            double * fidelity_full = iio_read_image_double_vec(
                        input_fidelity_string.c_str(),&width_f,&height_f,&d_f);
            for (int i=0; i<width_f*height_f*d_f; i++) {
                if (fidelity_full[i] != 0) {
                    fidelity_vec.push_back(i+1);
                    fidelity_vec.push_back(fidelity_full[i] -1);
                }
            }
        }
    
        int rows_fidelity = fidelity_vec.size()/2;
        double * fidelity = &fidelity_vec[0];

        /////////////Nystrom///////////////
        cout << "start Nystrom"<<endl;
        double * V_test = new double [N*m];
        double * D_test = new double [m];
        Nystrom N_obj;
        N_obj.Nystrom_algorithm(feature_vec,m,sigma,N,d,V_test,D_test);
        cout << "end Nystrom"<<endl;

        // output Nystrom to nystrom_result.txt
        ofstream nout;
        cout << "printing nystrom matrix to nystrom_result.txt\n";
        nout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/nystrom_result.txt");
        nout << setprecision(10) << fixed;
        if(nout.fail()) {
            cout << "But I couldn't open the file!!!\n";
            return 0;
        }
        for(int i=0; i<N; ++i) {
            for(int j=0; j<m; ++j) {
                nout << V_test[i*m + j] << " ";
            }
            nout << "\n";
        }

        /////////////MBO///////////////
        cout << "start semi-supervised MBO"<<endl;
        MBO_S mbo_object;
        double * index = new double[N];
        double * u0 = new double[N*n];
        // initialization
        mbo_object.MBO_algorithm(index, V_test, D_test, fidelity, u0,
                                 N, n, rows_fidelity, m,s, dt, mu,iterNum);
        cout << "end semi-supervised MBO"<<endl;
        // final index result
        double * final_index = new double [N+2];
        final_index[0]= N;
        final_index[1]= 1;
        for (int i = 0; i<N; i++) {
            final_index[i+2] = index[i];
        }
        //////////Output////////////////
        ofstream final;
        final.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt");
        final << final_index;
    }
    return 0;
}
