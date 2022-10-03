
#pragma once

#include <thread>
#ifndef UTILS_H
#define UTILS_H
void copyMatrix(double** origem, double** destino, int tamx, int tamy);
void flipMatrix(double** u, int tamx, int tamy);
void joinAll(std::thread* array, int tam);
#endif // UTILS_H
