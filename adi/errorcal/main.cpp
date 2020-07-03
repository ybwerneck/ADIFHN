#include <iostream>
#include "adifhn.h"
#include <sstream>
#include <string>

float calcError(FILE* referencia,FILE* solucao){


double vref;


while(fscanf(referencia, "%lf", &vref)==1){

        printf("%f \n",vref);




}

}
int main()

{
    FILE *referencia,*solucao;
    char endref[80];
    adifhn(2,-1,endref);

    return 0;

}
