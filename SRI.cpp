// SRI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread
#include "matrix.h"
#include <condition_variable>
#include "display.h"

Matrix* S, * R, * I;

//--------Parametros simulacao
int L = 4;
int t = 10;
float dt = 0.1, dx = 0.1;

int tam = L / dx;
//--------- Parametros modelo
float d11 = 0.01, d12 = 0.05, d22 = 0.01, d33 = 1;
float B = 1, A = 1;
float mi = 0.04, gama = 1, d = 1;
//---------

int n = t / dt;

void initValues() {

    S = new Matrix(tam, tam,0.99);
    I = new Matrix(tam, tam, 0.0);
    R = new Matrix(tam, tam, 0.0);
    

    I->operator()(tam/4 , tam / 4, 0.01);
    I->operator()(tam *1/3, tam * 1/ 3, 0.01);
    I->operator()(2*tam / 3, 2*tam / 3, 0.01);

}
float f(float S, float I) {

    return  - B * S * I;
}
float g(float S, float I) {

    return B * S * I  - (mi) * I;

}
float difussion(Matrix*ua,int x, int y,float coef) {

    //Difusão de Euler usada na aproximação de Ut+1 usada no calculo de Ru
    Matrix* u = ua;
    double difx, dify;
    float Y = (dt*coef) / (dx * dx);

    if (x != 0 && x != tam - 1)
        difx = (-2 * u->operator()(x, y) + u->operator()(x - 1, y) + u->operator()(x + 1, y)) * Y;
    else if (x == 0)
        difx = (-u->operator()(x, y) + u->operator()(x + 1, y)) * Y;
    else
        difx = (-u->operator()(x, y) + u->operator()(x - 1, y)) * Y;

    if (y != 0 && y != tam - 1)
        dify = (-2 * u->operator()(x, y) + u->operator()(x, y - 1) + u->operator()(x, y + 1)) * Y;
    else if (y == 0)
        dify = (-u->operator()(x, y) + u->operator()(x, y + 1)) * Y;
    else
        dify = (-u->operator()(x, y) + u->operator()(x, y - 1)) * Y;
    return difx + dify;

}
void step() {


    int i, j;
    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {

            float deltaS = dt * (f(S->operator()(i, j), I->operator()(i, j)) + difussion(S, i, j, d11) + difussion(I, i, j, d12));
            float deltaI = dt * (g(S->operator()(i, j), I->operator()(i, j)) + difussion(I, i, j, d22));
            float deltaR = dt * (mi * I->operator()(i, j) - d * R->operator()(i, j) + difussion(R, i, j, d33));
            S->operator()(i, j, deltaS + S->operator()(i, j));
            R->operator()(i, j, deltaR + R->operator()(i, j));
            I->operator()(i, j, deltaI + I->operator()(i, j));
        }
    }
}
void printMatrixtoFile(Matrix* I, Matrix* S, Matrix* R, FILE* Iarq, FILE* Rarq, FILE* Sarq) {

    for (int v = 0; v < tam; v++)
    {
        int dx = 1;
        for (int j = 0; j < tam; j++)
        {

            fprintf(Iarq, "%f ", I->operator()(v, j)); //att
            fprintf(Rarq, "%f ", R->operator()(v, j)); //att
            fprintf(Sarq, "%f ", S->operator()(v, j)); //att


        }
        fprintf(Iarq, "\n");
        fprintf(Sarq, "\n");
        fprintf(Rarq, "\n");
        
     };

    

    
}

int main()
{
      initValues();
    FILE* Iarq, * Sarq, * Rarq;
    char* filename = (char*)"I.txt";
    Iarq = fopen(filename, "w");
    char* filename2 = (char*)"S.txt";
    Sarq = fopen(filename2, "w");
    char* filename3 = (char*)"R.txt";
    Rarq = fopen(filename3, "w");

    for (int i = 0; i < n; i++) {

        step();
        printMatrixtoFile(I, S, R, Iarq, Rarq, Sarq);
       if(i!=n-1)
           {fprintf(Iarq, "\n\n");
           fprintf(Sarq, "\n\n");
           fprintf(Rarq, "\n\n");
       }

    }
    fclose(Iarq);
    fclose(Sarq);
    fclose(Rarq);
    saveGif((char*)"I.txt", (char*)"teste.gif", 0.1, 0.01);
    
}

