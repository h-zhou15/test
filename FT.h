#pragma once
#ifndef FT_H
#define FT_H
#include "pch.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
using namespace std;
typedef complex<double> complex_t;

int bitlen(int);
void DFT(vector<double> &x_n, vector<complex_t> &X_k);
void FFT(vector<complex_t> &x_n);
double format(double &c);
void IDFT(vector<complex_t> &X_k, vector<double> &x_n);
void IFFT(vector<complex_t> &x_n);
void resort(vector<complex_t> &x_n, int N);
int reverse_bit(int n, int len);
complex_t W(int k, int n, int N);

#endif