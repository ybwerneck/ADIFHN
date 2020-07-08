#include <iostream>
#include <sstream>
#include <string>
#include "adifhn.h"
#include "math.h"


double calcError(double a,double b) {

    double r = (a - b) / (1 - abs(a));
    return (pow(r,2));
}
double Error(char endref[80], char endsol[80]) {
    FILE* referencia, * solucao;
    referencia = fopen(endref, "r");  
    solucao = fopen(endsol, "r");
    double dxsol, dxref;
    double errot = 0;
    double ur,dxr,dyr,us,dxs,dys,dts,dtr;
    int k=1;
    int a;
    int i = 0;
    fscanf(referencia, "%lf", &dts);
    fscanf(solucao, "%lf", &dtr);
    if (abs(dts - dtr) > 0.001)
    {
        fclose(referencia);
        fclose(solucao);
        return 0;
    }    
    while (fscanf(solucao, "%lf %lf %lf", &us, &dxs, &dys)==3){
        do {  a = fscanf(referencia, "%lf %lf %lf", &ur, &dxr, &dyr); }
        while((dys != dyr || dxs!=dxr )&& a==3);
        if(a==3){
         errot+=calcError(us,ur);
		 i++;
        }
	
    }
    printf("%d pontos comparados\n", i);
    fclose(referencia);
    fclose(solucao);
    if (i < 100)
        return i;
    return sqrt(errot/i);
}
int main()

{
    int n = 1000;
    double dt = 0.3;
    char endref[80], endsol[80],endres[99];
	adifhn(n,dt, endref, true);
    system("cls");
    int zs[] =  {10,20,40,50,100,200};
   // sprintf(endref,"\resultados\\result_z1000.txt");
    double array[sizeof(zs)/sizeof(zs[0])];
    for (int i = 0; i < sizeof(zs) / sizeof(zs[0]); i++)
	{
		FILE* resultado;
		sprintf(endres, "resultadosErro\\resultado.txt");
		resultado = fopen(endres, "a");
        adifhn(zs[i], dt, endsol,true);
        array[i] = Error(endref, endsol);       
        fprintf(resultado,"%lf %lf\n", array[i], 1.0/zs[i]);
        fclose(resultado);
    }
    return 0;
    
}