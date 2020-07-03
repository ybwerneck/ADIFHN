#include <iostream>
#include <sstream>
#include <string>
#include "adifhn.h"
#include "math.h"


double calcError(double a,double b) {

    double r = (a - b) / (1 - abs(a));
    return sqrt(pow(r,2));
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
         errot+=calcError(us,ur);
        i++;
    }
    printf("%d pontos comparados\n", i);
    fclose(referencia);
    fclose(solucao);
    if (i < 100)
        return 0;
    return errot/i;
}
int main()

{
    int n =2000;
    double dt = 0.3;
    char endref[80], endsol[80];
    
    adifhn(n, dt, endref,true);
    system("cls");
    
    sprintf(endref,"result_z1000.txt");
    const int o = 100;
    double array[o-1];

    for (int i = 1; i < o; i++)
    {
        adifhn(i*10, dt, endsol,false);
        array[i-1] = Error(endref, endsol);
        printf("%f %d\n", array[i],i*10);
    }

    return 0;

}