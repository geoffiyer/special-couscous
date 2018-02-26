The Unupervised Graph-based Classification Algorithm 
********************************************************

*************
SUMMARY
*************

This program implements the unsupervised graph-based classification algorithm for image segmentation. Given an image, it computes the class labels of all the pixels and gives the segmentation of the image. It works for both hyperspectral and RGB images and it is parallelized so that it can work for video sequences too. 

This program is part of an IPOL publication:
http://www.ipol.im

*********
REFERENCE
*********

[1]  Zhaoyi Meng, Ekaterina Merkurjev, Alice Koniges and Andrea Bertozzi, “Hyperspectral Image Classification Using Graph Clustering Methods”, Image Processing Online, 2017.

***********
AUTHOR
***********

Zhaoyi Meng <mzhy@ucla.edu>
University of California, Los Angeles
Ekaterina Merkurjev <kmerkurev@math.msu.edu>
Michigan State University
Alice Koniges <aekoniges@lbl.gov>
Lawrence Berkeley National Laboratory
Andrea Bertozzi<bertozzi@math.ucla.edu>
University of California, Los Angeles


************
VERSION
************

Version 1, released on February 15, 2017

************
LICENSE
************

This program is free software: you can use, modify and/or redistribute it
under the terms of the simplified BSD License. You should have received a
copy of this license along this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

Copyright (C) 2015, Zhaoyi Meng <mzhy@ucla.edu>
All rights reserved.


******************
COMPILATION
******************

Required environment: Any unix-like system with a standard compilation 
environment (make and C and C++ compilers)

Required libraries: libtiff, cblas, clapack

Compilation instructions: run “make” to produce an executable “bin/b.exe”

*********
USAGE
*********

The program reads one input image, take some parameters and produce a classification model. The output is the segmented image. The meaning or the parameters is thoroughly discussed on the accompanying IPOL article. Usage instructions:

    <Usage>: ./bin/a.exe n  dt  mu  sigma  image  

    n           Number of classes 
    dt          The amount of time for which a diffusion takes place for
                 one iteration in the algorithm
    mu          The parameter before the fidelity term, larger lambda 
                 means you rely more on the known labels 
    sigma       The parameter in the definition of weights, larger sigma 
                 gives smaller weights between two nodes
    image       The RGB or hyperspectral image which you want to segment

Execution examples:

> ./bin/b.exe 4 0.15 0.15 0.1 ./data/cow.tiff 

******************
LIST OF FILES
******************
main.cpp             Main algorithm read the command line parameters and 
			run the algorithms
iio.c                Functions to read and write images
MBO_U.cpp            Functions of the unsupervised MBO algorithm
Nystrom.cpp          Functions of the Nystrom method. 
FeatureVec.cpp       Functions to compute the feature vectors for RGB images.


