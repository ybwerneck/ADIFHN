#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread
#include "matrix.h"
#include "adifhn.h"

double dx, dy;
double G;
double X;
double dt;
int Segundos;
int k;
double Y;
double L;
int tam;
double Cm;
double y, a, e;
bool fliped;
int n;
const int NT = 8;

Matrix initU()
{
	int i, j;
	Matrix u = Matrix(tam, tam, 0);
	for (i = 0; i < tam; i++)
	{
		for (j = 0; j < tam; j++)
		{
			u(i,j) = 0.0;

			int r = 0.2 / dx, c = (tam / 2);

			if (((i - c) * (i - c) + (j - c) * (j - c)) < r * r)
			{
				u(i,j) = 1;
			}
		}
	}
	return u;
}
double** iniciarMatrizRHS()
{
	double** matrix = (double**)malloc(sizeof(double*) * tam);
	int i, j;
	matrix[0] = (double*)calloc(sizeof(double), tam);
	matrix[1] = (double*)calloc(sizeof(double), tam);
	matrix[2] = (double*)calloc(sizeof(double), tam);

	for (i = 0; i < tam; i++)
	{
		matrix[0][i] = 1 + (2 * Y);
		matrix[1][i] = -(Y);
		matrix[2][i] = -(Y);
	}

	matrix[0][0] = 1 + Y;
	matrix[0][tam - 1] = 1 + Y;
	matrix[1][0] = 0.0;
	matrix[1][tam - 1] = -Y;
	matrix[2][0] = -Y;
	matrix[2][tam - 1] = 0.0;

	return matrix;
}
double* getVetorFx(Matrix &Ua, Matrix &g, int x)
{
	double* u = (double*)calloc(sizeof(double), tam);
	int i=0;
	for (i = 0; i < tam; i++) {
		
		double U = Ua(x, i),UXM1,UXP1;
		UXM1 = (x != 0) ? Ua(x-1,i) : 0, UXP1 = (x != tam - 1) ? Ua(x+1,i) : 0;
		
		int y = i;
		double Ru = dt * ((U * (1 - U) * (U - a) - g(x,y)) / e);
		if (x != 0 && x != tam - 1)
			u[i] = U * (1 - 2 * Y) + Y * UXP1 + Y * UXM1 + Ru;
		else if (x == 0)
			u[i] = U * (1 - Y) + Y * UXP1 + Ru;

		else u[i] = U * (1 - Y) + Y * UXM1 + Ru;
	}
	
	return u;
}
void updateG(Matrix &u, Matrix &g, int n, int bgx, int tamx, int tamy)
{
	int i, j;

	for (i = bgx; i < tamx; i += n)
	{
		for (j = 0; j < tamy; j++)
		{
			g(i,j, g(i,j)+ dt * 0.5 * (u(i,j) - (g(i,j) * y)));
		}
	}
}

void updateGt(Matrix &u, Matrix &g, int tam)
{
	int i;
	const int n = NT;
	std::thread array[n];
	for (int i = 0; i < n; i++)
	{
		array[i] = std::thread(updateG, std::ref(u), std::ref(g), n, i, tam, tam);
	}
	joinAll(array, n);
}
double* solve_tridiagonal(double x[], const size_t N, const double a[], const  double b[], const  double c[])
{
	size_t n;

	double* const cprime = (double*)malloc(sizeof(double) * N);

	cprime[0] = c[0] / b[0];
	x[0] = x[0] / b[0];

	for (n = 1; n < N; n++)
	{
		double m = 1.0f / (b[n] - a[n] * cprime[n - 1]);
		cprime[n] = c[n] * m;
		x[n] = (x[n] - a[n] * x[n - 1]) * m;
	}

	for (n = N - 1; n--> 0.0; )
		x[n] = x[n] - cprime[n] * x[n + 1];

	free(cprime);
	return x;
}
void print(int n, int k, double kt, bool report) {
	if (report) {
		int o = kt == -1 ? k : (kt) / dt;
		system("cls");
		printf("DT: %.10f DX:%.10f TAM= %d \n Quadro %d \ %d   completed %f", dt, dx, tam, n, o, ((double)n / (o)) * 100);
	}
}

void calc(Matrix &u, Matrix &ua, int beg, int n, int tam, double** rhs, Matrix &g) {
	for (int i = beg; i < tam; i += n)
	{
		double* Fx=getVetorFx(ua, g, i);
		Fx = solve_tridiagonal(Fx, tam, rhs[1], rhs[0], rhs[2]);
		
		if (ua.isFliped())
			u.setLine(i,Fx);
		else
			u.setColumn(i,Fx);
		
	}
}

Matrix step (Matrix &ua, double** rhs, Matrix  &g) {
	Matrix u=Matrix(tam,tam,0);
	if (ua.isFliped())
		u.flip();
	const int n = NT;
	std::thread array[n];

	for (int i = 0; i < n; i++)
	{
		array[i] = std::thread(calc, std::ref(u), std::ref(ua), i, n, tam, rhs, std::ref(g));
	}
	joinAll(array, n);
	updateGt(u, g, tam);
	
	return u;
}
void adifhn(double dtP, double dxP, double dtAlvo, char c[], bool report)
{

	FILE* arquivo;
	char* filename = c;
	arquivo = fopen(filename, "w");
	dx = dxP, dy = dxP;
	L = 1;
	tam = L / dx;
	G = 25;
	X = 100;
	dt = dtP;
	Segundos = 2;
	k = Segundos / dt;
	Y = (G / X) * ((dt) / (2 * dx * dx));
	Cm = 0.0;
	y = 0.5, a = 0.1, e = 0.01;
	Matrix u= initU();
	Matrix g= Matrix(tam,tam,0);
	double** rhs = iniciarMatrizRHS();
	fliped = false;
	int n;
	for (n = 0; n < k; n++)
	{
		std::thread printer = std::thread(print, n, k, dtAlvo, report);
		printer.detach();
		if (dtAlvo == -1 || dtAlvo <= n * dt) {
			if (dtAlvo != -1)
				fprintf(arquivo, "%f \n", dt * n);

			for (int v = 0; v < tam; v++)
			{
				for (int j = 0; j < tam; j++)
				{
					double x = v * dx;
					double y = j * dx;
					double t = n * dt;

					if (dtAlvo == -1)
						fprintf(arquivo, "%.6f %.6f %.6f %f \n", t, x, y, u(v, j));
					else
						fprintf(arquivo, "%f  %f  %f\n", u(v, j), dx * v, dy * j);
				}

				fprintf(arquivo, "\n");
			}

			int i;
			if (dtAlvo == -1)
			{
				fprintf(arquivo, "\n\n");
			}
		}
		if (dtAlvo <= n * dt && dtAlvo != -1)
			break;

		u = step(u, rhs, g);
	//	u.print();

		u.flip();
		g.flip();
		u = step(u, rhs, g);
	//	u.print();
		u.flip();
		g.flip();
	}

	
	fclose(arquivo);

	printf("%s foi salvo!", filename);
}