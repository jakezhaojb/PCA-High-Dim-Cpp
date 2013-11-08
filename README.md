PCA_High_Dim
============

>**Designer:** Junbo Zhao, Wuhan University, Working in Tsinghua National lab of intelligent images and documents processing.   
**Email:** zhaojunbo1992chasing@gmail.com	      +86-18672365683

**Introduction:**     
This package includes a small and convenient lib for Matrix operations, and a Principle 
Component Analysis (PCA) program. This PCA project is designed to be suitable to High 
Dimensional situations, using some linear algebra tricks to transform the covariance matrix and 
achieving a faster and less memory-cost way for PCA.

**Platform:**     
This program is tested on VS2010, 32bit, Win7 system. I cannot guarantee it could be adopted on 
other platforms.    
For Windows users, you can just add the .cpp and .h files in your project and compile it.    
For Linux users, I will soon update this program with a makefile.   
For Mac users, I am sorry that this program is and will not be tested on Mac system. If you can 
successfully compile this program on Mac, please feel free to contact me！ 

**Files:**     
In detail, "eigen.cpp" as well as "eigen.h" combined finishes computing eigenvalues and corresponding eigenvectors;     
"matrix.cpp" and "matrix.h" in charge of all the matrix operations;      
"pca.cpp" implements PCA.

**Usage:**      
1. You should prepare your .txt file of the feature matrix before implementation.    
The specific rules of .txt file can be found in the sample.txt.     
2. View the main function (in pca.cpp) to see how to exploit PCA.
