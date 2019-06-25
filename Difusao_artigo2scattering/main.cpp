
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;
double L;
void print_vector(double *v, int j)
{
    cout<<endl;
    for(int i=0; i<j; i++)
    {
        cout << v[i] << endl;
    }
    cout<<endl;
}

double R (double u)
{
    return u*(1-u)*(u-0.8);
//collison -0.8
//scattering -0.1
}
double CondicaoInicial(double x)
{
    if((x>=L/3)&& x<=2*(L/3))
    {
        return 1;
    }
    else
    {
    return 0;
    }
}
int main()
{
    FILE *arquivo;
    double d, dx, dt, alfa, aux;
    int n, j;
    L=10.0;
    arquivo= fopen("arquivos/colision.txt","w");
    if(arquivo)
    {
        d=0.001;
        dx=0.04;
        dt=0.05;
        cout<<"dx= "<<dx<<endl;
        cout<<"dt= "<<dt<<endl;
        j= (int)(L/dx); //números de pontos em x
        n= 4000; //números de pontos em t
        cout<<"j: "<<j<<" n: "<<n<<endl;
        alfa= (d*dt)/pow(dx,2);
        double u_novo[j], u_antigo[j];



        for(int i=0; i<=j; i++)
        {

            //define condição inicial e imprime os valores pra t=0
            u_antigo[i]=CondicaoInicial(i*dx);
        }

        for(int t=0; t<=n; t++)//for no tempo
        {

            // cout<<"Tempo : "<<t<<endl;
            for (int x=0; x<j-1; x++) //for no espaço
            {
                u_novo[x]= ((R(u_antigo[x])*dt) + u_antigo[x] + (alfa*(u_antigo[x+1]-(2*u_antigo[x])+u_antigo[x-1])));

            }
            u_novo[0]= u_novo[1];
            u_novo[j-1]= u_novo[j-2];


            for(int i=0; i<=j; i++)
            {
                u_antigo[i]=u_novo[i];
            }
            if(t%100==0)
                for (int x=0; x<=j; x++) //for no espaço
                {

                    //  cout<<"Espaço : "<<x*dx<<endl;

                    fprintf(arquivo, "%f ", u_novo[x]);

                    fprintf(arquivo, "%f ", (dx*x));
                    fprintf(arquivo, "%d ", t);
                    fprintf(arquivo, "\n");

                }
        }
        fclose(arquivo);
    }
    return 0;
}

