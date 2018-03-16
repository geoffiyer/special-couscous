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
#include<cassert>
#include "MBO_S.h"
#include "Nystrom.h"
#include <string>
#include <vector>
#include "FeatureVec.h"
#include <cblas.h>
extern "C" {
#include "iio.h"
}

// #define PRINT_NYSTROM

using namespace std;

// this function is for writting the result
// the first two elements of A are the dimentions
std::ostream& operator<<(std::ostream& out, double *A);

int main(int argc, char * argv[]){

    //////////////////////////////////////////////
    // I prefer to have parameters hardcoded in
    /////////////////////////////////////////////
    // Parameters to adjust are all here
    int dim1 = 3;          // dim of first modality
    int dim2 = 1;          // dim of second modality
    double dt;       // MBO algorithm step size
    double mu;       // fidelity parameter. Set through file input
    bool twoNorm = 1;      // Norm is either 2norm or cos angle thing
                           // corresponds to true and false, resp
    bool nystrom = 0;      // True means do the nystrom. False means load from file

    // These only are used in the case that we do nystrom via c++ (instead of load)
    double sigma = 0.8;    // weight matrix scaling parameter
    double lambda = 1;     // scaling between the two modalities
    // (I think high # favors modality number 2)
    
    int m;             // Nystrom size
    int n;             // num classes

    // Work with whatever matlab just put out.
    string input_image_string = "/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/dataFromMatlab.tiff";
    string input_fidelity_string = "/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/fidelityFromMatlab.tiff";
    string input_evec_string("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/");
    string input_evec_string_angle("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab_Angle/"); 

    // Umbrella data inputs
    // n = 6;             // num classes
    // string input_image_string = "/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/umbrella_both.tiff";
    // string input_fidelity_string = "/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/fidelity_umbrella_6class.tiff";
    // string input_evec_string("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/Umbrella_matlab/");
    // string input_evec_string_angle("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/Umbrella_angle/"); 
   
    // // DFC data inputs
    // n = 6;             // num classes    
    // string input_image_string("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/DFC_both.tiff");
    // string input_fidelity_string("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/fidelity_DFC_6classSMOOTH.tiff");
    // string input_evec_string("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/DFC_matlab/");
    // string input_evec_string_angle("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/DFC_angle/");

    // // Jadeplant data input
    // n = 8;
    // string input_image_string = "/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/jadeplant_both_nonlocal.tiff";
    // string input_fidelity_string = "/home/gsiyer/Schoolwork/Chanussot/MBO_Code/data/fidelity_jadeplant_8class.tiff";
    // string input_evec_string("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/Jadeplant_matlab/");
    // string input_evec_string_angle("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Input/Jadeplant_angleNonlocal/");

    /////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////
    // Leftover from Gloria code (inputs come through terminal here)
    /////////////////////////////////////////////////////////////////
    // int n = atoi(argv[1]);
    // double dt = atof(argv[2]);
    // double mu = atof(argv[3])*100.0;
    // double sigma = atof(argv[4]);
    // string input_image_string(argv[5]);
    // string input_fidelity_string;
    // int dim1 = atoi(argv[6]);    // dimension of first modality
    // int dim2 = atoi(argv[7]);    // dimension of second modality
    // double lambda = atof(argv[8]);
    // input_fidelity_string = argv[9];
    // bool twoNorm = true;
    /////////////////////////////////////////////////////////////////
    
    double iterNum= 500;
    int s = 3;  // number of diffusions steps before each thresholding.
    int windowSize = 2;

    if(nystrom)
        m = 100;    
    else {
        ifstream nmin;
        string temp;
        if(twoNorm)
            temp = input_evec_string + "nm.txt";
        else
            temp = input_evec_string_angle + "nm.txt";
        nmin.open(temp.c_str());
        assert(!nmin.fail());
        nmin >> n;
        nmin >> m;
        nmin >> mu;
        nmin >> dt;
    }

    /////////////Input image///////////////
    // int width;
    // int height;
    // int d;
    // double * input_image;
    // input_image = iio_read_image_double_vec(input_image_string.c_str(),
    //     &width,&height,&d);
    // assert( dim1+dim2 == d );
    // cout << "Width = " << width
    //      << "\nHeight = " << height
    //      << "\ndim1 = " << dim1
    //      << "\ndim2 = " << dim2 << "\n";
        
    // // add a small number to all the pixel values in case there are a lot of 0s
    // for (int i=0; i<width*height*d; i++) {
    //     input_image[i] += 0.1;
    // }

    // ofstream outImage;
    // cout << "Loaded image is printed to image.txt\n";
    // outImage.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/imageHeader.txt");
    // outImage << width*height << "\n" << dim1 << "\n" << dim2;
    // outImage.close();
    // outImage.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/image.txt");
    // for(int i=0; i<width*height; i++) {
    //     for(int j=0; j<d; j++) {
    //         outImage << input_image[i*d + j] << " ";
    //     }
    //     outImage << "\n";
    // }

    // double * feature_vec;
        
    // if (d<4) {
    //     // for RGB image, rescale some parameters, so that it is
    //     // more consistent with the hyperspectral images
    //     sigma = sigma*100.0;
    //     mu = mu*10.0;
    //     // for RGB image, need to read frame by frame and build the feature vector
    //     double * input_image_split = new double [width*height*d];
            
    //     for (int i=0; i<width*height; i++) {
    //         for (int l=0; l<d; l++) {
    //             input_image_split[width*height*l+i] = input_image[d*i+l];
    //         }
    //     }
            
    //     FeatureVec Feature_obj;
            
    //     double * kernel = new double[(2*windowSize+1)*(2*windowSize+1)];
    //     for (int i=0; i<(2*windowSize+1)*(2*windowSize+1); i++) {
    //         kernel[i] = 0;
    //     }
    //     Feature_obj.makeKernel(kernel,windowSize);
            
    //     int N1p = height + 2 * windowSize;
    //     int N2p = width + 2 * windowSize;
    //     double* padimage = new double[N1p*N2p*3];
    //     for (int i=0; i<N1p*N2p*3; i++) {
    //         padimage[i] = 0;
    //     }
    //     Feature_obj.padarray(input_image_split,padimage,
    //         windowSize,height,width);
            
    //     int size_feature = width*height*(2*windowSize+1)*(2*windowSize+1)*3;
    //     feature_vec = new double[size_feature];
            
    //     Feature_obj.makeFeatureVec(feature_vec,padimage,kernel,
    //         windowSize,height,width);
            
    //     d = (2*windowSize+1)*(2*windowSize+1)*3;
    //     delete [] padimage;
    // } else {
    //     cout << "Hyperspectral Image"<<endl;
    //     feature_vec = input_image;
    // }
    ////////////// End Input Image ///////////
    
    // int N = width*height;  // without image one does not need width height
    
    // input fidelity
    ofstream outImage;
    int width_f;
    int height_f;
    int d_f;
    vector<double> fidelity_vec;

    double * fidelity_full = iio_read_image_double_vec(
        input_fidelity_string.c_str(),&width_f,&height_f,&d_f);
    outImage.close();
    outImage.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt");
    cout << "Fidelity Stats:\nwidth = " << width_f << "\nheight: " << height_f
         << "\ndim: " << d_f << "\n"
         << "Printing fidelity loaded into fidelity.txt\n";
    cout << "Fidelity const: " << mu << "\n";
    cout << "dt const: " << dt << "\n";
            
    for (int i=0; i<width_f*height_f*d_f; i++) {
        outImage << fidelity_full[i] << "\n";
        if (fidelity_full[i] != 0) {
            fidelity_vec.push_back(i);
            fidelity_vec.push_back(fidelity_full[i]-1);
        }
    }
    // assert(width_f == width && height_f == height);
    int width = width_f;
    int height = height_f;
    int N = width*height;
    
    outImage.close();
    
    int rows_fidelity = fidelity_vec.size()/2;
    double * fidelity = &fidelity_vec[0];

    ////////////////
    //  I was nervous about initializing fidelity like above, because if fidelity_vec
    //  dies we have a hanging pointer. Then I realized this will never happen.
    //  They go out of scope at the same time. I think.....
    //  Anyways, the code below is less scary to me, but makes a copy.
    ////////////////
    // double * fidelity = new double [fidelity_vec.size()];
    // for(int i=0; i < rows_fidelity*2; ++i) {
    //     fidelity[i] = fidelity_vec[i];
    // }
    ////////////////

    ////// If you want to print the fidelity read to console /////////////
    // cout << "rows_fidelity: " << rows_fidelity << "\n";
    // for(int i=0; i<fidelity_vec.size(); ++i) {
    //     cout << fidelity_vec[i] << "\n";
    // }
    ///////////////////////////////////////////////////////
    
    /////////////Nystrom///////////////
    double * V_test = new double [N*m];
    double * D_test = new double [m];

    if(nystrom) {
        // m = 100;
        // Nystrom N_obj;
        // cout << "start Nystrom"<<endl;
        // N_obj.Nystrom_algorithm(feature_vec,m,sigma,N,dim1,dim2,lambda,
        //     V_test,D_test,twoNorm);
        // cout << "end Nystrom"<<endl;
        // // output Nystrom to nystrom_result.txt
        // ofstream nout;
        // cout << "printing nystrom matrix to nystrom_V.txt and nystrom_D.txt\n";
        // nout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/nystrom_V.txt");
        // nout << setprecision(13);
        // if(nout.fail()) {
        //     cout << "But I couldn't open the file!!!\n";
        //     return 0;
        // }
        // for(int i=0; i<N; ++i) {
        //     for(int j=0; j<m; ++j) {
        //         nout << V_test[i*m + j] << " ";
        //     }
        //     nout << "\n";
        // }
        // nout.close();
        // nout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/nystrom_D.txt");
        // for(int i=0; i<m; ++i) {
        //     nout << D_test[i] << "\n";
        // }
    }
    else {
        ifstream in;
        cout << "loading evecs from Input/nystrom_V.txt\n";
        string temp;
        if(twoNorm)
            temp = input_evec_string + "nystrom_V.txt";
        else
            temp = input_evec_string_angle + "nystrom_V.txt";
        in.open(temp.c_str());
        assert(!in.fail());
        for(int i=0; i<N; ++i) {
            for(int j=0; j<m; ++j) {
                in >> V_test[i*m + j];
            }
        }
        in.close();
        if(twoNorm)
            temp = input_evec_string + "nystrom_D.txt";
        else
            temp = input_evec_string_angle + "nystrom_D.txt";
        in.open(temp.c_str());
        for(int i=0; i<m; ++i) {
            in >> D_test[i];
        }
        cout << "evecs loaded\n";

#ifdef PRINT_NYSTROM        
        // output Nystrom to nystrom_result.txt
        ofstream nout;
        cout << "printing nystrom matrix to nystrom_V.txt and nystrom_D.txt\n";
        nout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/nystrom_V.txt");
        assert(!nout.fail());
        nout << setprecision(13);
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
        nout.close();
        nout.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/nystrom_D.txt");
        assert(!nout.fail());
        for(int i=0; i<m; ++i) {
            nout << D_test[i] << "\n";
        }
#endif

    }
        
    /////////////MBO///////////////
    cout << "start semi-supervised MBO"<<endl;
    MBO_S mbo_object;
    double * index = new double[N];
    double * u0 = new double[N*n];
    // initialization
    mbo_object.MBO_algorithm(index, V_test, D_test, fidelity, u0,
        N, n, rows_fidelity, m, s, dt, mu,iterNum);
    cout << "end semi-supervised MBO"<<endl;
    // final index result
    double * final_index = new double [N+2];
    final_index[0]= N;
    final_index[1]= 1;
    for (int i = 0; i<N; ++i) {
        final_index[i+2] = index[i];
    }
    //////////Output////////////////
    ofstream final;
    final.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt");
    assert(final.is_open());
    final << final_index;
    final.close();

    return 0;
}


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
