/*
  ģ����Բ��ͶӰֵ���ɹ���
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

//���ڿ��Բ�ȡ����ͶӰֵ
//double** Projection_IO(string path){}
//

//ģ��ƽ��ͶӰ������¾�����Բ��ͶӰֵ
void simulation_p(double **pValue, double Attenuation, double a, double b) {
	//	double DisInterval = 1;//ƽ����ֱ�߼��
	double LongAxis = max(a, b);//�ж���Բ���ᣬ����ɨ��ķ�Χ

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
	int theta_len = sizeof(pValue[0]) / sizeof(double); //��Ӧtheta����ֵ
	int distance_len = len / theta_len;//��ӦDistance����ֵ

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

//���ø���Ҷ��Ƭ����ͶӰֵת��Ϊ��ά��ʽ��Ƶ��ͶӰֵ
void  FourierSlice(double **pValue)
{
	vector<vector<complex_t>> X_k = ConvertRS(pValue);
	for (int Theta = 0; Theta < 180; Theta++)
	{
		
	}
}