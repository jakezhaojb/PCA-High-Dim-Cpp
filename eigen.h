#ifndef _EIGEN_H_
#define _EIGEN_H_

#include <valarray>
#include "matrix.h"
#include <math.h>
#include <iostream>

using namespace std;

int QR(valarray<double>& main_diag, valarray<double>& minor_diag, Matrix &q, double accuracy, int l);
int Householder(Matrix m, Matrix &q, valarray<double>& main_diag, valarray<double>& minor_diag);

#endif
