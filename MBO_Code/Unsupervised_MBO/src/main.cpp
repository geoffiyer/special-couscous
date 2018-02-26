// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017 Zhaoyi Meng <mzhy@ucla.edu>
// All rights reserved.

#include <iomanip>
#include <time.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <assert.h>
#include "Nystrom.h"
#include "MBO_U.h"
#include "FeatureVec.h"
extern "C" {
#include "iio.h"
}
using namespace std;

void normalize(double * A, int rows, int cols){
    for(int i = 0; i<rows; i++){
        double sum(0);
        for(int j = 0; j<cols; j++){
            sum = sum + A[i*cols+j]*A[i*cols+j];
        }
        double norm(0);
        norm = sqrt(sum);
        for(int j = 0; j<cols; j++){
            A[i*cols+j] = A[i*cols+j]/norm;
        }
    }
}

// the first two elements of A are the dimentions
std::ostream& operator<<(std::ostream& out, double *A){
    int rows;
    int cols;
    rows = A[0];
    cols = A[1];
    
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++)
            out<<A[i*cols+j+2]<<"\t";
        out<<"\r\n";
    }
    out<<"\r\n\n";
    
    return out;
}


void printEigenfunctions(std::ostream& out, const double *V,
    int N, int m) {
    
    for(int i=0; i<N; i++) {
        for(int j=0; j<m; j++) {
            out << V[i*m + j] << "\t";
        }
        out << "\r\n";
    }

    out << "\r\n\n";
        
}
int main(int argc, const char * argv[]) {
    if (argc<6) {
        cout << "Not enough input parameters"<<endl;
        printf("Usage: %s n c1 dt sigma input_image_string \n", argv[0]);
    } else {
        srand(time(NULL));
        
        int n = atoi(argv[1]);        //number of classes
        double dt = atof(argv[2]);
        double lambda = atof(argv[3])*100000.0;
        double sigma = atof(argv[4]);
        string input_img_string(argv[5]);  //input image
        int windowSize = 2;
        
        // read in image of tiff format
        int width;
        int height;
        int d;
        double * input_image;
        input_image = iio_read_image_double_vec(input_img_string.c_str(),
            &width,&height,&d);

        cout << "Width = " << width
             << "\nHeight = " << height
             << "\nd = " << d << "\n";

        int M = min(100,width*height);
        int m = min(100,width*height);
        
        double * feature_vec;
        if (d<4) {
            // for RGB image, rescale some parameters so that it is more
            // consistant with the hyperspectral image
            sigma = sigma * 100.0;
            lambda = lambda / 100.0;
            // for RGB image, need to read frame by frame
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
            Feature_obj.padarray(input_image_split,padimage,windowSize,
                height,width);
            
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
        
        // Nystrom
        cout << "start Nystrom algorithm"<<endl;
        double * V = new double [N*m];
        double * D = new double [m];
        
        Nystrom N_obj;
        N_obj.Nystrom_algorithm(feature_vec,m,sigma,N,d,V,D);
        cout << "end Nystrom algorithm"<<endl;

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
                if(V[i*m+j]>0)
                    nout << " ";
                nout << V[i*m + j] << " ";
            }
            nout << "\n";
        }

        // Mumford_Shah
        cout << "start unsupervised MBO algorithm"<<endl;
        // TODO: why do we normalize here?
        normalize(feature_vec, N, d);
        
        // the first 2 element of final_index are the dimensions
        double *final_index = new double[N+2];
        
        MBO_U MS_obj;
        MS_obj.Mumford_Shah(final_index,feature_vec,V,D,N,n,M,d,dt,lambda);
        cout << "end unsupervised MBO algorithm"<<endl;
        ofstream final;
        final.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt");
        assert(final.is_open());
        final << final_index;
        final.close();
        
        if (d<4) {
            delete [] feature_vec;
        }
        delete [] input_image;
        delete [] V;
        delete [] final_index;
    }
    
    return 0;
}
