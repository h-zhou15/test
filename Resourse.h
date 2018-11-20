#pragma once
#ifdef RESOURSE_H
#define RESOURSE_H

#include"pch.h"
#include"Projection.h"

template<class T>
void VectorIndexSwap2D(const T & _x)
{	
	T _y = _x;
	_x.clear();
	for (int i = 0; i < _x.size(); i++)
	{
		for (int j = 0; j < _x[0].size(); j++)
			_x[j][i] = _y[i][j];
	}
}

#endif // RESOURSE_H
