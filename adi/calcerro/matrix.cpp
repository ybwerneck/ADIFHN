#include "matrix.h"
using namespace std;

// Constructor for Any Matrix
Matrix::Matrix(unsigned rowSize, unsigned colSize, double initial) {
	m_rowSize = rowSize;
	m_colSize = colSize;
	m_matrix.resize(rowSize);
	for (unsigned i = 0; i < m_matrix.size(); i++)
	{
		m_matrix[i].resize(colSize, initial);
	}
}

Matrix::Matrix(const Matrix& B)
{
	this->m_colSize = B.getCols();
	this->m_rowSize = B.getRows();
	this->m_matrix = B.m_matrix;

}

Matrix::~Matrix() {

}


// Returns value of given location when asked in the form A(x,y)
double& Matrix::operator()(const unsigned& rowNo, const unsigned& colNo)
{

	if (!fliped)
		return this->m_matrix[rowNo][colNo];
	else
		return this->m_matrix[colNo][rowNo];


}
void Matrix::operator()(const unsigned& rowNo, const unsigned& colNo, double Value)
{
	if(!fliped)
		this->m_matrix[rowNo][colNo]=Value;
	else 
		this->m_matrix[colNo][rowNo] = Value;

}


// No brainer - returns row #
unsigned Matrix::getRows() const
{
	return this->m_rowSize;
}

// returns col #
unsigned Matrix::getCols() const
{
	return this->m_colSize;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose()
{
	Matrix Transpose(m_colSize, m_rowSize, 0.0);
	for (unsigned i = 0; i < m_colSize; i++)
	{
		for (unsigned j = 0; j < m_rowSize; j++) {
			Transpose(i, j) = this->m_matrix[j][i];
		}
	}
	return Transpose;
}

// Prints the matrix beautifully
void Matrix::print() const
{
	cout << "Matrix: " << endl;
	for (unsigned i = 0; i < m_rowSize; i++) {
		for (unsigned j = 0; j < m_colSize; j++) {
			printf("[%.2f]",m_matrix[i][j]);
		}
		cout << endl;
	}
}

bool Matrix::flip()
{
	fliped = !fliped;
	return fliped;
}

bool Matrix::isFliped()
{
	return fliped;
}

void Matrix::setColumn(const int col,double* val) 
{
	for (unsigned i = 0; i < m_rowSize; i++) {
		m_matrix[col][i] = val[i];
	}
}

void Matrix::setLine(const int lin, double* val)
{
	for (unsigned i = 0; i < m_rowSize; i++) {
		m_matrix[i][lin] = val[i];
	}
}

