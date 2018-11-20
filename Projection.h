
#pragma once
#ifndef PROJECTION_H
#define PROJECTION_H
#include "pch.h"
/*
  模拟椭圆的投影值生成过程
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>


typedef std::complex<double> complex_t;

using namespace std;
//模拟平行投影束情况下均匀椭圆的投影值
void simulation_p(double **pValue, double Attenuation, double a, double b);

vector<vector<complex_t>>  FourierSlice_fft(double **pValue);
vector<vector<complex_t>> FourierSlice_dft(double **pValue);
vector<vector<complex_t>> ConvertRS_fft(double **pValue);	 //为fft产生匹配序列
vector<vector<double>> ConvertRS_dft(double **pValue);	 //为dft产生匹配序列

void Filter(vector<vector<complex_t>> &F_w);

vector<vector<complex_t>> BackProjection(vector<vector<complex_t>> &P_F);

//double** Projection_IO(string path);//待添加

template<class T>
void VectorIndexSwap2D(T & _x)
{
	T _y = _x;
	_x.clear();
	for (int i = 0; i < _x.size(); i++)
	{
		for (int j = 0; j < _x[0].size(); j++)
			_x[j][i] = _y[i][j];
	}
}

#endif