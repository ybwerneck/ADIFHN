#pragma once
//
//  matrix.h
//
//  Created by Furkanicus on 12/04/15.
//  Copyright (c) 2015 Furkan. All rights reserved.
// Learnrope.com

#ifndef __EE_242_Project_2__matrix__
#define __EE_242_Project_2__matrix__

#include <stdio.h>
#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

using std::vector;
using std::tuple;
class Matrix {
private:
	unsigned m_rowSize;
	unsigned m_colSize;
	vector<vector<double> > m_matrix;
	bool fliped=false;
public:
	Matrix(unsigned, unsigned, double);
	Matrix(const char*);
	Matrix(const Matrix&);
	~Matrix();
	Matrix transpose();

	double& operator()(const unsigned&, const unsigned&);
	void operator()(const unsigned&, const unsigned&, double);
	void print() const;
	bool flip();
	bool isFliped();
	void setColumn(const int, double*);
	unsigned getRows() const;
	unsigned getCols() const;
	void setLine(const int col, double* val);

};
#endif /* defined(__EE_242_Project_2__matrix__) */
