
#pragma once

#include <thread>
#ifndef UTILS_H
#define UTILS_H
void flipmatrixt(double** u, int tamx, int tamy);
void flipmatrixa(double** u, int tamx, int tamy);
void joinAll(std::thread* array, int tam);
#endif // UTILS_H
