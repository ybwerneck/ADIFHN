#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread

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

double** initU()
{
	int i, j;
	double** U = (double**)calloc(sizeof(double*), tam);
	for (i = 0; i < tam; i++)
	{
		U[i] = (double*)calloc(sizeof(double), tam);
		for (j = 0; j < tam; j++)
		{
			U[i][j] = 0.0;

			int r = 0.2 / dx, c = (tam / 2);

			if (((i - c) * (i - c) + (j - c) * (j - c)) < r * r)
			{
				U[i][j] = 1;
			}
		}
	}
	return U;
}

double** initG()
{
	int i, j;
	double** G = (double**)calloc(sizeof(double*), tam);
	for (i = 0; i < tam; i++)
	{
		G[i] = (double*)calloc(sizeof(double), tam);
		for (j = 0; j < tam; j++)
		{
			G[i][j] = 0;
		}
	}
	return G;
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
double* getVetorFx(double** Ua, double** g, int x)
{
	double* u = (double*)calloc(sizeof(double), tam);
	int i;
	for (i = 0; i < tam; i++) {
		double U = Ua[x][i], UXM1 = (x != 0) ? Ua[x - 1][i] : 0, UXP1 = (x != tam - 1) ? Ua[x + 1][i] : 8;
		int y = i;
		double Ru = dt*((U * (1 - U) * (U - a) - g[x][y]) / e);
		if (x != 0 && x != tam - 1)
			u[i] = U * (1 - 2 * Y) + Y * UXP1 + Y * UXM1 + Ru;
		else if (x == 0)
			u[i] = U * (1 - Y) + Y * UXP1 + Ru;

		else u[i] = U * (1 - Y) + Y * UXM1 + Ru;
	}
	
	return u;
}
void updateG(double** g, double** V, int n, int bgx, int tamx, int tamy)
{
	int i, j;

	for (i = bgx; i < tamx; i += n)
	{
		for (j = 0; j < tamy; j++)
		{
			g[i][j] = g[i][j] + dt * 0.5 * (V[i][j] - (g[i][j] * y));
		}
	}
}

void updateGt(double** V, double** g, int tam)
{
	int i, j;
	const int n = NT;
	std::thread array[n];
	for (int i = 0; i < n; i++)
	{
		array[i] = std::thread(updateG, g, V, n, i, tam, tam);
	}
	joinAll(array, n);
}
void solve_tridiagonal(double x[], const size_t N, const double a[], const  double b[], const  double c[])
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

	for (n = N - 1; n-- > 0.0; )
		x[n] = x[n] - cprime[n] * x[n + 1];

	free(cprime);
}
void print(int n, int k, double kt, bool report) {
	if (report) {
		int o = kt == -1 ? k : (kt) / dt;
		system("cls");
		printf("DT: %.10f DX:%.10f TAM= %d \n Quadro %d \ %d   completed %f", dt, dx, tam, n, o, ((double)n / (o)) * 100);
	}
}

void calc(double** u, double** ua, int beg, int n, int tam, double** rhs, double** g) {
	for (int i = beg; i < tam; i += n)
	{
		u[i] = getVetorFx(ua, g, i);
		solve_tridiagonal(u[i], tam, rhs[1], rhs[0], rhs[2]);
	}
}

double** step(double** ua, double** rhs, double** g) {
	double** u;
	u = (double**)malloc(sizeof(double*) * tam);

	const int n = NT;
	std::thread array[n];
	for (int i = 0; i < n; i++)
	{
		array[i] = std::thread(calc, u, ua, i, n, tam, rhs, g);
	}
	joinAll(array, n);
	updateGt(u, g, tam);

	for (int i = 0; i < tam; i++)
	{
		delete[] ua[i];
	}

	free(ua);
	return u;
}
void adifhn(double dtP, double dxP, double dtAlvo, char c[], bool report)
{
	FILE* arquivo;

	char* filename = c;
	arquivo = fopen(filename, "w");

	dx = dxP, dy = dxP;
	L = 1.1;

	tam = L / dx;

	G = 25;
	X = 100;
	dt = dtP;
	Segundos = 1;
	k = Segundos / dt;
	Y = (G / X) * ((dt) / (2 * dx * dx));
	Cm = 0.0;
	y = 0.5, a = 0.1, e = 0.01;
	double** u;
	double** g;
	double** rhs = iniciarMatrizRHS();
	u = initU();
	g = initG();
	fliped = false;
	int n;
	double** L = initU();
	
	for (n = 0; n < k; n++)
	{
		std::thread printer = std::thread(print, n, k, dtAlvo, report);
		printer.detach();

		u = step(u, rhs, g);
		flipmatrixa(u, tam, tam);
		flipmatrixa(g, tam, tam);
		fliped = true;

		u = step(u, rhs, g);
		flipmatrixa(g, tam, tam);
		flipmatrixa(u, tam, tam);
		fliped = false;

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
						fprintf(arquivo, "%.6f %.6f %.6f %f \n", t, x, y, u[v][j]);
					else
						fprintf(arquivo, "%f  %f  %f\n", u[v][j], dx * v, dy * j);
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
	}

	free(u);
	free(g);

	fclose(arquivo);

	printf("%s foi salvo!", filename);
}