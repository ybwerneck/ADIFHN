#include <iostream>
#include <sstream>
#include <string>
#include "adifhn.h"
#include "math.h"
#include "display.h"


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
    if (abs(dts - dtr) > 0.0001)
    {
        fclose(referencia);
        fclose(solucao);
        return -dtr;
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
        return -i;
    return sqrt(errot/i);
}
int main()

{
    double dtT = -1;
    char endref[80], endsol[80], endres[99];                             
    sprintf(endref, "resultados\\resultref.txt");
    sprintf(endres, "resultadosErro\\resultado.txt");
    double dx = 0.1, dt = dx * dx;
     adifhn(dt, dx, dtT, endref, true);
    system("cls");
    exibirGif(endref,dx, dt);
    return 0;
    FILE* resultado;
	resultado=fopen(endres, "w");
	fclose(resultado);

    double zs[] = {0.1,0.05,0.02,0.01,0.005,0.002,0.001};
    double array[sizeof(zs) / sizeof(zs[0])];
    for (int i = 0; i < sizeof(zs) / sizeof(zs[0]); i++)
    {
        double dx = zs[i], dt = dx * dx;
        sprintf(endsol, "resultados\\result_DX=%f.txt", dx);
        resultado = fopen(endres, "a");
        adifhn(dt, dx, dtT, endsol, true);
        array[i] = Error(endref, endsol);
        fprintf(resultado, "%lf %lf\n", array[i], dx);
        char resultadofoto[200];
		sprintf(resultadofoto, "resultados\\result_DX=%f.png", dx);
        saveFoto(endsol,resultadofoto, dx, dt);
        fclose(resultado);
        
    }
	return 0;

}