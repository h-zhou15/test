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
#include <complex>

#define PI 3.14159265354

//���ڿ��Բ�ȡ����ͶӰֵ
//double** Projection_Read(string path){}
//

//ģ��ƽ��ͶӰ������¾�����Բ��ͶӰֵ
void simulation_p(double **pValue, double Attenuation, double a, double b) {

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

vector<vector<complex_t>> ConvertRS_fft(double **pValue)
{
	int len = sizeof(pValue) / sizeof(double);
	int theta_len = sizeof(pValue[0]) / sizeof(double); //��Ӧtheta����ֵ
	int distance_len = len / theta_len;//��ӦDistance����ֵ

	vector<vector<complex_t>> P_k;
	
	for (int i = 0; i < theta_len; i++)
	{
		for (int j = 0; j < distance_len; j++)
		{
			P_k[i].push_back(complex_t(pValue[j][i], 0));
		}
	}
	return P_k;			//Ϊp(theta,distance)��ʽ												
}

vector<vector<double>> ConvertRS_dft(double **pValue)
{
	int len = sizeof(pValue) / sizeof(double);
	int theta_len = sizeof(pValue[0]) / sizeof(double); //��Ӧtheta����ֵ
	int distance_len = len / theta_len;//��ӦDistance����ֵ

	vector<vector<double>> X_n_dft;

	for (int i = 0; i < theta_len; i++)
	{
		for (int j = 0; j < distance_len; j++)
		{
			X_n_dft[i].push_back(pValue[j][i]);
		}
	}
	return X_n_dft;					//Ϊp��theta,distance)��ʽ
}

//���ø���Ҷ��Ƭ����ͶӰֵת��Ϊ��ά��ʽ��Ƶ��ͶӰֵ
vector<vector<complex_t>>  FourierSlice_fft(double **pValue)
{
	vector<vector<complex_t>> FST_X_k = ConvertRS_fft(pValue);

	for (int Theta = 0; Theta < FST_X_k.size(); Theta++)
	{
		FFT(FST_X_k[Theta]);			 //����Ҷ��Ƭ����
	}

	return FST_X_k;
}

vector<vector<complex_t>> FourierSlice_dft(double **pValue)
{
	vector<vector<double>> x_n = ConvertRS_dft(pValue);

	vector<vector<complex_t>> FST_X_k;

	for (int Theta = 0; Theta < x_n.size(); Theta++)
	{
		DFT(x_n[Theta], FST_X_k[Theta]);			//FST_X_k(theta,w)����ʽ��wΪ��ɢ����µ�kֵ
	}
	return FST_X_k;		   //δ�˲�ǰ��ͶӰֵ
}

//�˲������� 
void  Filter(vector<vector<complex_t>> &X_k)
{
	//ѡ�ü򵥵��˲�����	   |w|
	for (int theta = 0; theta < X_k.size(); theta++)
	{
		for (int k = 0; k < X_k[theta].size(); k++)
		{
			X_k[theta][k] = complex_t(X_k[theta][k].real()*k, X_k[theta][k].imag()*k); //�˲�
		}
	}
	//��������Ҷ���任���õ��˲����ͶӰ
	for (int theta = 0; theta < X_k.size(); theta++)
	{
		IFFT(X_k[theta]);	//�˲����ͶӰֵ
	}
}

vector<vector<complex_t>> BackProjection(vector<vector<complex_t>>& P_F)
{
	VectorIndexSwap2D<vector<vector<complex_t>>>(P_F);	//����ΪP_f(w,theta)�ĸ�ʽ
	vector<vector<complex_t>> _fxy;

	//��ͶӰ
	for (int x = 0; x < P_F.size(); x++)
	{
		for (int y = 0, theta = 0; theta < P_F[theta].size(); y++, theta++)
		{
			_fxy[x][y] += P_F[x][theta];	  
		}
	}
	return _fxy;
}
