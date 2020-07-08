#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"
#define n 8

void joinAll(std::thread* array, int tam) {

    for (int i = 0; i < tam; i++)
    {
        array[i].join();

    }
}
void dump_vectore(double* u, int t)
{
    int j;

    printf("[");
    for (j = 0; j < t; j++)


        printf("(i %d j% d - %f \n ", t, j, u[j]);

    printf("]");


}
void dump_vector(double* u, int tam)
{
    int j;
    printf("[");

    for (j = 0; j < tam; j++)
        printf(" %f ", u[j]);


    printf("]");
}
void dump_matrix(double** u, int tamx, int tamy)
{
    int i, j;
    for (i = 0; i < tamx; i++)
    {
        printf("[");
        for (j = 0; j < tamy; j++)
            printf(" %f ", u[i][j]);
        printf("] \n");
    }
}

void copyMatrixe(double** origem, double** destino, int bgx, int tamx, int tamy)
{
    int i, j;

    for (i = bgx; i < tamx; i += n)
    {
        for (j = 0; j < tamy; j++)
        {

            destino[i][j] = origem[i][j];

        }
    }

}


void copyMatrix(double** origem, double** destino, int tamx, int tamy)
{
	const int nt = n;
	std::thread array[nt];
    int i;
	for (i = 0; i < n; i++)
	{
		array[i] = std::thread(copyMatrixe, origem, destino, i, tamx, tamy);
	}
	joinAll(array, nt);

}

void allocmatrix(double** u, double** L, int k, int tamx) {
    int i, j; for (i = k; i < tamx; i += n)
    {
        L[i] = (double*)calloc(sizeof(double), tamx);

    }
}
void unallocmatrix(double** u, double** L, int k, int tamx) {
    int i, j; for (i = k; i < tamx; i += n/2)
    {
        delete[] L[i];

    }
}
void flipmatrix(double** u, double** L, int k, int tamx, int tamy) {
    int i, j; for (i = k; i < tamx; i += n)
    {
        for (j = k; j < tamy; j++)
        {
            L[i][j] = u[j][i];
        }
    }
}
void flipmatrixt(double** u, int tamx, int tamy)
{

    double** L = (double**)calloc(sizeof(double*), tamy);
    int i, j;
    std::thread array[n];
    for (int i = 0; i < n; i++)
    {
        array[i] = std::thread(allocmatrix, u, L, i,tamx);
    }
    joinAll(array, n);

    for (int i = 0; i < n; i++)
    {
        array[i] = std::thread(flipmatrix, u, L, i, tamx, tamy);
    }
    joinAll(array, n);

    for (i = 0; i < n; i++)
    {
        array[i] = std::thread(copyMatrixe, L, u, i, tamx, tamy);
    }
    joinAll(array, n);

    for (i = 0; i < tamx; i++)
    {

        delete[] L[i];
    }
    delete[] L;
}



void flipmatrixa(double** u, int tamx, int tamy)
{
    const int nt = n;
    std::thread array[nt];
    int i, j;
    double** L = (double**)calloc(sizeof(double*), tamx);
    for (i = 0; i < tamx; i++)
    {
        L[i] = (double*)calloc(sizeof(double), tamx);
        for (j = 0; j < tamy; j++)
        {
            L[i][j] = u[j][i];
        }
    }
    
    for (i = 0; i < n; i++)
    {
         array[i] = std::thread(copyMatrixe, L,u, i, tamx, tamy);
    }
    joinAll(array, nt);
    
    for (i = 0; i < tamx; i++)
    {

        delete[] L[i];
    }
    delete[] L;
}
