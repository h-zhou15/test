
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

typedef std::complex<double> complex_t;
using namespace std;

//模拟平行投影束情况下均匀椭圆的投影值
void simulation_p(double **pValue, double Attenuation, double a, double b);
void  FourierSlice(double **pValue);
vector<vector<complex_t>> ConvertRS(double **pValue);
//double** Projection_IO(string path);//待添加

#endif