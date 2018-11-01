
#pragma once
#ifndef PROJECTION_H
#define PROJECTION_H
#include "pch.h"
/*
  ģ����Բ��ͶӰֵ���ɹ���
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

typedef std::complex<double> complex_t;
using namespace std;

//ģ��ƽ��ͶӰ������¾�����Բ��ͶӰֵ
void simulation_p(double **pValue, double Attenuation, double a, double b);
void  FourierSlice(double **pValue);
vector<vector<complex_t>> ConvertRS(double **pValue);
//double** Projection_IO(string path);//�����

#endif