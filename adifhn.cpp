/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread
#include "matrix.h"
#include "adifhn.h"
#include <mutex>
#include <condition_variable>
#include <atomic>         // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT
#include <chrono>
using namespace std::chrono_literals;
#include <chrono>
using namespace std::chrono_literals;

double dx, dy;
double G;
int npar=0;

double X;
double dt;
int Segundos;
int k;
double Y;
double L;
int tam;
double Cm;
double y, a, e;
int o = 0;
const int NT = 10;
bool* flaguG, * flagCal, * flaguGkill, * flagCalkill;
Matrix* u, * ua, * g;
std::mutex mutex;
std::condition_variable cUg, cCal,cs,cs2;
int done(0);
Matrix* initU()
{
	int i, j;
	Matrix* u = new Matrix(tam, tam, 0);
	for (i = 0; i < tam; i++)
	{
		for (j = 0; j < tam; j++)
		{
			u->operator()(i, j, 0.0);

			int r = (L *0.2) / dx, c = (tam / 2);

			if (((i - c) * (i - c) + (j - c) * (j - c)) < r * r)
			{
				u->operator()(i, j, 1);
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
double* getVetorFx(int x)
{
	double* Fx;
	Fx = (double*)malloc(sizeof(double) * tam);

	int i = 0;
	for (i = 0; i < tam; i++) {
	
		double U = ua->operator()(x, i), UXM1, UXP1;
		UXM1 = (x != 0) ? ua->operator()(x - 1, i) : 0, UXP1 = (x != tam - 1) ? ua->operator()(x + 1, i) : 0;
		// que beleza de sintaxe em c++
		int y = i;
		double Ru = dt* ((U * (1 - U) * (U - a) -g->operator()(x, y)) / e);
		if (x != 0 && x != tam - 1)
			Fx[i] = U * (1 - 2 * Y) + Y * UXP1 + Y * UXM1 + Ru;
		else if (x == 0)
			Fx[i] = U * (1 - Y) + Y * UXP1 + Ru;

		else Fx[i] = U * (1 - Y) + Y * UXM1 + Ru;
	}

	return Fx;
}

double* solve_tridiagonal(double x[], const size_t N, const double a[], const  double b[], const  double c[])
{
	size_t n;
	double* cprime;

	cprime = (double*)malloc(sizeof(double) * tam);

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

	delete cprime;
	return x;
}
void print(int k, double kt, bool report) {
	while (report){
		int n = npar;
		int o = kt == -1 ? k : (kt) / dt;
		system("cls");
			printf("\nDT: %.10f DX:%.10f TAM= %d \n Quadro %d \ %d   completed %f", dt, dx, tam, n, o, ((double)n / (o)) * 100);
			std::this_thread::sleep_for(1000ms);
	}
}

void calc(int beg, int tam, double** rhs) {
	while (flagCalkill[beg]==false) {

		std::unique_lock<std::mutex> lock(mutex); // aquire mutex lock as required by condition variable
		cCal.wait(lock, [beg] {return (flagCal[beg]||flagCalkill[beg]); }); // thread will suspend here and release the lock if the expression does not return true
		lock.unlock();
		if (flagCalkill[beg] == true)
			break;
		for (int i = beg; i < tam; i += NT)
		{
			double* Fx = getVetorFx(i);
			Fx = solve_tridiagonal(Fx, tam, rhs[1], rhs[0], rhs[2]);
			
			for (int j = 0; j < tam; j++)
			{
				u->operator()(i, j, Fx[j]);
				g->operator()(i, j, g->operator()(i, j) + dt * 0.5 * (u->operator()(i, j) - (g->operator()(i, j) * y)));
			}
			delete[] Fx;

		}

		lock.lock();
		done++;

		lock.unlock();

		flagCal[beg] = false;
		cs.notify_all();

	}
}

void calctridiag() {
	done = 0;

	for (int k = 0; k < NT; k++)
	{
		flagCal[k] = true; 

	}
	cCal.notify_all();
	std::unique_lock<std::mutex> lock(mutex); // aquire mutex lock as required by condition variable
	cs.wait_for(lock,100ms, [] {return(done >= NT); }); // thread will suspend here and release the lock if the expression does not return true
}
void adifhn(double dtP, double dxP, double dtAlvo, char c[], bool report)
{
	FILE* arquivo;
	char* filename = c;
	arquivo = fopen(filename, "w");
	
	//iniciação de variaveis
	dx = dxP, dy = dxP;
	L = 1.1;
	tam = L / dx;
	G = 25;
	X = 100;
	dt = dtP;
	Segundos = 2;
	k = Segundos / dt;
	Y = (G / X) * ((dt) / (2 * dx * dx));
	Cm = 0.0;
	y = 0.5, a = 0.1, e = 0.01;
	g = new Matrix(tam, tam, 0);
	u = initU();
	ua = initU(); 
	double** rhs = iniciarMatrizRHS();
	
	//iniciação das threads 
	std::thread threadsCalc[NT];
	flagCal = (bool*)malloc(sizeof(bool) * NT);
	flagCalkill = (bool*)malloc(sizeof(bool) * NT);
		for (int i = 0; i < NT; i++)
	{
		flagCal[i] = false;
	}

		for (int i = 0; i < NT; i++)
	{
		flagCalkill[i] = false;
	}

		for (int i = 0; i < NT; i++)
	{
		threadsCalc[i] = std::thread(calc, i, tam, rhs);


	}
	std::thread printer = std::thread(print, k, dtAlvo, report);
	printer.detach();



	
	
	int n;
	for (n = 0; n < k; n++)  //LOOP PRINCIPAL
	{
		npar = n;
		if ((n%1==0)&(dtAlvo == -1 || dtAlvo <= n * dt)) {
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
						fprintf(arquivo, "%.6f %.6f %.6f %f %f \n", t, x, y, ua->operator()(v, j), g->operator()(v, j));
					else
						fprintf(arquivo, "%f  %f  %f\n", ua->operator()(v, j), dx * v, dy * j);
				}

				fprintf(arquivo, "\n");
			}

			int i;
			if (dtAlvo == -1)
			{
				fprintf(arquivo, "\n\n");
			}
		}
	    if (dtAlvo <= n * dt && dtAlvo != -1) {
			
			break;

		}
		calctridiag(); //40 ms por loop 1100x 1100
		u->flip();
		ua->flip();			// CalculaFx e depois o usa para Calcular a tridiag, ulitliza NT threads
		g->flip();			//Flip -> inversão logica no index das matrixes 
		std::swap(u, ua);
		calctridiag();
		u->flip();
		ua->flip();
		g->flip();
		std::swap(u, ua);

	}
	for (int i = 0; i < NT; i++)
	{
		flagCalkill[i] = true;
		cCal.notify_all();
		threadsCalc[i].join();

	}

	fclose(arquivo);

	printf("%d", o);
}

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread
#include "matrix.h"
#include "adifhn.h"
#include <mutex>
#include <condition_variable>
#include <atomic>         // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT
#include <chrono>
using namespace std::chrono_literals;
#include <chrono>

double dx, dy;
double G;
int npar = 0;

double X;
double dt;
int Segundos;
int k;
double Y;
double L;
int tam;
double Cm;
double y, a, e;
int o = 0;
const int NT = 6;
bool* flaguG, * flagCal, * flaguGkill, * flagCalkill;
Matrix* u, * ua, * g;
std::mutex mutex, mutex2, mutex3;
std::condition_variable cUg, cCal, cs, cs2;
std::atomic<int> done;
std::atomic<bool> report;
Matrix* initU()
{
	int i, j;
	Matrix* u = new Matrix(tam, tam, 0);
	for (i = 0; i < tam; i++)
	{
		for (j = 0; j < tam; j++)
		{
			u->operator()(i, j, 0.0);

			int r = (L * 0.2) / dx, c = (tam / 2);

			if (((i - c) * (i - c) + (j - c) * (j - c)) < r * r)
			{
				u->operator()(i, j, 1);
			}
		}
	}
	return u;
}

bool areDone(bool* flag, int t) {
	for (int i = 0; i < t; i++)
		if (flag[i])
			return false;
	return true;
}
void print(int k, double kt) {
	while (report) {
		int n = npar;
		int o = kt == -1 ? k : (kt) / dt;
		system("cls");
		printf("\nDT: %.10f DX:%.10f TAM= %d \n Quadro %d \ %d   completed %f", dt, dx, tam, n, o, ((double)n / (o)) * 100);
		std::this_thread::sleep_for(300ms);
	}
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
double difusaoEuler(int x, int y) {

	//Difusão de Euler usada na aproximação de Ut+1 usada no calculo de Ru
	Matrix* u = ua;
	double difx, dify;
	if (x != 0 && x != tam - 1)
		difx = (2 * u->operator()(x, y) - u->operator()(x - 1, y) - u->operator()(x + 1, y)) * Y;
	else if (x == 0)
		difx = (u->operator()(x, y) - u->operator()(x + 1, y)) * Y;
	else
		difx = (u->operator()(x, y) - u->operator()(x - 1, y)) * Y;
	if (y != 0 && y != tam - 1)
		dify = (2 * u->operator()(x, y) - u->operator()(x, y - 1) - u->operator()(x, y + 1)) * Y;
	else if (y == 0)
		dify = (u->operator()(x, y) - u->operator()(x, y + 1)) * Y;
	else
		dify = (u->operator()(x, y) - u->operator()(x, y - 1))* Y;
	return difx + dify;

}
double R(double U, double g) {
	return  dt * ((U * (1 - U) * (U - a) - g) / e);
}
double* getVetorFx(int x)
{
	double* Fx;
	Fx = (double*)malloc(sizeof(double) * tam);

	int y = 0;
	for (y = 0; y < tam; y++) {

		double U = ua->operator()(x, y), UXM1, UXP1,G=g->operator()(x,y),Ru;
		UXM1 = (x != 0) ? ua->operator()(x - 1, y) : 0, UXP1 = (x != tam - 1) ? ua->operator()(x + 1, y) : 0;
		if (true) {   // true para ultilizar uma aproximação de U ulitlizando o método de Euler no calculo de Ru

			U = (R(U,G) + difusaoEuler(x, y)) * dt + U;
		}
		Ru = R(U,G);
		if (x != 0 && x != tam - 1)
			Fx[y] = U * (1 - 2 * Y) + Y * UXP1 + Y * UXM1 + Ru;
		else if (x == 0)
			Fx[y] = U * (1 - Y) + Y * UXP1 + Ru;

		else Fx[y] = U * (1 - Y) + Y * UXM1 + Ru;
	}

	return Fx;
}
double* solve_tridiagonal(double x[], const size_t N, const double a[], const  double b[], const  double c[])
{
	size_t n;
	double* cprime;

	cprime = (double*)calloc(sizeof(double), tam);

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

	delete[] cprime;
	return x;
}
void calc(int beg, int tam, double** rhs) {
	while (flagCalkill[beg] == false) {
		cs.notify_one();
		std::unique_lock<std::mutex> lock(mutex); 
		//Espera hora de iniciar novo passo
		do cCal.wait_for(lock,10ms, [beg] {return (flagCal[beg] || flagCalkill[beg]); }); 
		while (!(flagCal[beg] || flagCalkill[beg]));
		lock.unlock();
		if (flagCalkill[beg] == true)
			break; //Finaliza thread com flag
		
		//Passo
		for (int i = beg; i < tam; i += NT)
		{
			double* Fx = getVetorFx(i);
			Fx = solve_tridiagonal(Fx, tam, rhs[1], rhs[0], rhs[2]);

			for (int j = 0; j < tam; j++)
			{
				u->operator()(i, j, Fx[j]);
				g->operator()(i, j, g->operator()(i, j) + dt * 0.5 * (u->operator()(i, j) - (g->operator()(i, j) * y)));
			}
			delete[] Fx;

		}
		//-----------Passo
		flagCal[beg] = false; //Avisa que acabou a execução do passo. 
		cs.notify_one();



	}
}
void calctridiag() {
	int k;

	//Seta a flag que indica para as threads que devem executar um passo 	for (done=0,k = 0; k < NT; k++)
	for(k=0;k<NT;k++){
		flagCal[k] = true;


	}
	cCal.notify_all();
    std::unique_lock<std::mutex> lock(mutex); 
	//Espera a execução terminar sem busywait
	do cs.wait_for(lock, 100ms, [] {return(areDone(flagCal, NT));}); 
		while (!areDone(flagCal,NT));
}
void adifhn(double dtP, double dxP, double dtAlvo, char c[], bool printB)
{
	//Inicialização de variaveis
	FILE* arquivo;
	char* filename = c;
	arquivo = fopen(filename, "w");
	dx = dxP, dy = dxP;
	L = 1.1;
	tam = L / dx;
	G = 25,	X = 100;
	dt = dtP;
	Segundos = 2;
	k = Segundos / dt;
	Y = (G / X) * ((dt) / (2 * dx * dx));
	Cm = 0.0,y = 0.5, a = 0.1, e = 0.01;
	double** rhs = iniciarMatrizRHS();
	g = new Matrix(tam, tam, 0);
	u = initU();
	ua = initU();
	flagCal = (bool*)malloc(sizeof(bool) * NT);
	flagCalkill = (bool*)malloc(sizeof(bool) * NT);
	//-------------------Inicialização de variaveis


	//Inicialização de Threds
	std::thread threadsCalc[NT];
	for (int i = 0; i < NT; i++)
	{
		flagCal[i] = false;
	}
	for (int i = 0; i < NT; i++)
	{
		flagCalkill[i] = false;
	}
	for (int i = 0; i < NT; i++)
	{
		threadsCalc[i] = std::thread(calc, i, tam, rhs);


	}
	report = printB;
	std::thread printer = std::thread(print, k, dtAlvo);
	//----------------Inicialização de Threds

	int n;
	for (n = 0; n < k; n++)
	{
		npar = n;
		//Printa matrix no arquivo
		if ((n % 1 == 0) & (dtAlvo == -1 || dtAlvo <= n * dt)) { 
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
						fprintf(arquivo, "%.6f %.6f %.6f %f %f \n", t, x, y, ua->operator()(v, j), g->operator()(v, j));
					else
						fprintf(arquivo, "%f  %f  %f\n", ua->operator()(v, j), dx * v, dy * j);
				}

				fprintf(arquivo, "\n");
			}

			int i;
			if (dtAlvo == -1)
			{
				fprintf(arquivo, "\n\n");
			}
		}
		//Para a execução apos atingir o tempo desejado
		if (dtAlvo <= n * dt && dtAlvo != -1) {break;}
		
		//Loop principal


		//t->t+1/2 
		calctridiag(); // Calcula FX e resolve a matrix tridiagonal, ultiliza ua com referencia e assigna valores obtidos em U
		u->flip();		//Inverte a matrix em uma inversão lógica 
		ua->flip();
		g->flip();
		std::swap(u, ua); //Passa os Valores de U para a variavel de referencia

		//t+1/2->t+1  com a matrix invertida
		calctridiag();
		u->flip();
		ua->flip();
		g->flip();
		std::swap(u, ua);
		//----Loop principal

	}
	//Encerrar Threads e fechar arquivo
	for (int i = 0; i < NT; i++)
	{
		flagCalkill[i] = true;
		cCal.notify_all();
		threadsCalc[i].join();

	}

	fclose(arquivo);
	report = false;
	printer.join();
	printf("%d", o);
}