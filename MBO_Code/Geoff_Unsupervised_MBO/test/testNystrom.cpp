#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<algorithm>
#include<cmath>

extern "C" {
#include<cblas.h>
#include"iio.h"
}
using namespace std;

extern "C"{
    double ddot_( const int *, const double *, const int *, const double *,
                 const int * );
    double dnrm2_(const int *, const double *, const int *);
}

double scale_stdev(double* A, double* B, int n1, int n2);
    
int main() {
    
    std::ofstream out;
    out.open("results.txt");
    if(out.fail()) {
        std::cout << "Ack open fail\n";
        return 0;
    }
    
    double* super_test;
    int width;
    int height;
    int d;
    super_test = iio_read_image_double_vec("../data/test.tiff",&width,&height,&d);
    int dim1 = d;
    int dim2 = 0;
    double sigma = 1;

    int N = width*height;
    int m = N/2;
    cout << "Width = " << width
         << "\nHeight = " << height
         << "\ndim1 = " << dim1
         << "\ndim2 = " << dim2 << "\n";
    
    std::vector<int> myvector;
    for (int i=0; i<N; ++i) myvector.push_back(i);
//    std::random_shuffle ( myvector.begin(), myvector.end() );
    
    // creating Wxx (m by m) and Wxy (m by N-m)    
    double* Wxx2 = new double [m*m];
    double* Wxx = new double[m*m];
    double* differenceVec = new double[d];
    int inca = 1, incb = 1;
    for(int i = 0; i<m; i++){
        int temp_1 = myvector[i]*d;
        for(int j = 0; j<m; j++){
            int temp_2 = myvector[j]*d;
            for(int k=0; k<d; k++) {
                differenceVec[k] = super_test[temp_1+k] - super_test[temp_2+k];
            }
            Wxx [i*m+j] = dnrm2_(&dim1, differenceVec, &inca);
            Wxx2[i*m+j] = dnrm2_(&dim2, differenceVec+dim1, &incb);
        }
    }
    
    // Wxy , m by N-m
    double* Wxy  = new double [m*(N-m)];
    double* Wxy2 = new double [m*(N-m)];
    for(int j = 0; j<N-m; j++){
        int temp_1 = myvector[j+m]*d;
        for(int i = 0; i<m; i++){
            int temp_2 = myvector[i]*d;
            for(int k=0; k<d; k++) {
                differenceVec[k] = super_test[temp_1+k] - super_test[temp_2+k];
            }
            Wxy [i*(N-m)+j] = dnrm2_(&dim1, differenceVec, &inca);
            Wxy2[i*(N-m)+j] = dnrm2_(&dim2, differenceVec+dim1, &incb);
        }
    }
    
    scale_stdev(Wxx , Wxy , m*m, m*(N-m));
    scale_stdev(Wxx2, Wxy2, m*m, m*(N-m));

    for(int i=0; i<m*m; i++)
        Wxx[i] = std::max(Wxx[i], Wxx2[i]);

    for(int i=0; i<m*(N-m); i++)
        Wxy[i] = std::max(Wxy[i], Wxy2[i]);

    delete[] Wxx2;
    delete[] Wxy2;
    delete[] differenceVec;

    for(int i=0; i<m*m; i++)
        Wxx[i] = exp(-Wxx[i]/sigma);
    for(int i=0; i<m*(N-m); i++)
        Wxy[i] = exp(-Wxy[i]/sigma);
    
    for(int i=0; i<m; i++) {
        for(int j=0; j<N; j++) {
            if(j < m)
                out << Wxx[i*m+j] << " ";
            else
                out << Wxy[i*(N-m) + (j-m)] << " ";
        }
        out << "\n";
    }

    return 0;
}

double scale_stdev(double* A, double* B, int n1, int n2) {
    double mean = 0;

    for(int i=0; i<n1; i++)
        mean += A[i];
    for(int i=0; i<n2; i++)
        mean += B[i];

    mean = mean/(n1+n2);
    double stdev = 0;

    for(int i=0; i<n1; i++)
        stdev += (A[i] - mean)*(A[i] - mean);
    for(int i=0; i<n2; i++)
        stdev += (B[i] - mean)*(B[i] - mean);
    stdev = sqrt(stdev/(n1+n2));

    if(stdev < 1e-11)
        return stdev;

    for(int i=0; i<n1; i++) 
        A[i] /= stdev;
    for(int i=0; i<n2; i++)
        B[i] /= stdev;

    return stdev;
}
