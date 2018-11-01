/*
  模拟椭圆的投影值生成过程
*/
#include "pch.h"
#include "Projection.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "FT.h"
#define PI 3.14159265354

typedef std::complex<double> complex_t;

using namespace std;

//后期可以采取读入投影值
//double** Projection_IO(string path){}
//

//模拟平行投影束情况下均匀椭圆的投影值
void simulation_p(double **pValue, double Attenuation, double a, double b) {
	//	double DisInterval = 1;//平行束直线间隔
	double LongAxis = max(a, b);//判断椭圆长轴，决定扫描的范围

	for (int Theta = 0; Theta < 180; Theta++) {
		for (int Distance = -(LongAxis+0.5); Distance <= int(LongAxis+0.5); Distance += 1) {
			double r2 = pow(a*cos(PI*Theta / 180), 2) + pow(b*sin(PI*Theta / 180), 2);
			if (r2 - Distance * Distance >= 0) {
				double PQ = 2 * a*b*sqrt(r2 - Distance * Distance) / r2;
				pValue[Distance + (int)(LongAxis+0.5)][Theta] = PQ * Attenuation;
			}
			else
				pValue[Distance + (int)(LongAxis + 0.5)][Theta] = 0;
		}
	}
}

vector<vector<complex_t>> ConvertRS(double **pValue)
{
	int len = sizeof(pValue) / sizeof(double);
	int theta_len = sizeof(pValue[0]) / sizeof(double); //对应theta索引值
	int distance_len = len / theta_len;//对应Distance索引值

	vector<vector<complex_t>> P_k;
	
	for (int i = 0; i < distance_len; i++)
	{
		for (int j = 0; j < theta_len; j++)
		{
			P_k[i].push_back(complex_t(pValue[i][j], 0));
		}
	}
	return P_k;

}

//利用傅里叶切片定理将投影值转换为二维形式的频域投影值
void  FourierSlice(double **pValue)
{
	vector<vector<complex_t>> X_k = ConvertRS(pValue);
	for (int Theta = 0; Theta < 180; Theta++)
	{
		
	}
}