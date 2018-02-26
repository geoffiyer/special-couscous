
#include<iostream>
#include<fstream>
#include<cmath>
extern "C"{
    double dnrm2_(const int *, const double *, const int *);
}

using namespace std;

int main() {

    int N;         // num pixels
    int dim1;      // dim of 1st modality
    int dim2;      // dim of 2nd modality
    int d;         // total dim 
    double* im;    // 1D array representing image
    double* K;     // 1D array representing classification

    double scaling1;  // stdev scaling of 1st modality
    double scaling2;  // stdev scaling of 2nd modality
    
    ifstream in;
    in.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/imageHeader.txt");
    if(in.fail()) {
        cout << "Egad there is an error opening imageHeader\n";
        return 0;
    }
    in >> N >> dim1 >> dim2;
    d = dim1 + dim2;
    in.close();

    in.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/scaling.txt");
    if(in.fail()) {
        cout << "Egad there is an error opening scaling\n";
        return 0;
    }
    in >> scaling1 >> scaling2;
    in.close();

    im = new double[N*d];
    in.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/image.txt");
    if(in.fail()) {
        cout << "Egad there is an error opening image\n";
        return 0;
    }
    for(int i=0; i<N*d; ++i)
        in >> im[i];
    in.close();

    K = new double[N];
    in.open("/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt");
    if(in.fail()) {
        cout << "Egad there is an error opening classification result\n";
        return 0;
    }
    int temp;
    in >> temp;
    in >> temp;
    for(int i=0; i<N; ++i) 
        in >> K[i];
    in.close();
    
    double error = 0;
    double* differenceVec = new double[d];
    int inca = 1, incb = 1;
    for(int i = 0; i<N; i++){
        int temp_1 = i*d;
        for(int j = i+1; j<N; j++){
            int temp_2 = j*d;
            if(K[i+2] != K[j+2]) {
                for(int k=0; k<d; k++) {
                    differenceVec[k] = im[temp_1+k] - im[temp_2+k];
                }
                error += max(dnrm2_(&dim1, differenceVec, &inca)/(scaling1),
                    dnrm2_(&dim2, differenceVec+dim1, &incb)/(scaling2))/N;
            }
        }
    }
    delete[] differenceVec;

    cout << "Error: " << error << "\n";

    delete[] im;
    delete[] K;
    
}
