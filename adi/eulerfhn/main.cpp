#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;


//20
//10
//13
//17
//+
//36
const double dx=0.1,dy=0.1;
const double L=3;
const int nx=L/dx + 2;
const double Ly=3;
const int ny=Ly/dy + 2;

double Ge=25,X=100,Cm=1.0;
double dt=0.001;
double e=0.01, a=0.1;
double G = (Ge/X)*((dt)/(dx*dx));

int t=2;
int n=t/dt;


double f(double U,double g)
{

    return (U*(1 - U)*(U - a) - g)/e;

}
double difusao(double Va[nx][ny],double Ka[nx][ny],int x,int y)
{
    return  (((G)*Va[x-1][y]-2*(G)*Va[x][y]+(G)*Va[x+1][y])) +   (((G)*Va[x][y-1]-2*(G)*Va[x][y]+(G)*Va[x][y+1]));
}
double g(double v, double k)
{
    double gk=1; //veruufcar
    double y=0.5;
    return k+(dt*(v-(y*k)));
}
int main()
{
    printf("%f",G);


    FILE *arquivo;
    FILE *arquivo2;

    arquivo= fopen("result.txt","w");

    arquivo2= fopen("result2.txt","w");

    double V[nx][ny];
    double K[nx][ny];
    double Va[nx][ny];
    double Ka[nx][ny];


    //condição inicial
    for (int i=0; i<=nx-1; i++)
    {
        for (int j=0; j<=ny-1; j++)
        {
            Va[i][j] = 0;
            Ka[i][j] =0;

            int r=nx/10,c=(nx/2)-1;
            if(((i-c)*(i-c) + (j-c)*(j-c)) < r*r )
            {
                Va[i][j]=1;

            }



        }
    }
    for (int k=0; k<=n; k++) //for no tempo
    {


        //loop fora da borda
        for(int x=1; x<nx-1; x++)
        {
            for(int y=1; y<nx-1; y++)
            {

                V[x][y]= dt*2*(f(Va[x][y],Ka[x][y]))+difusao(Va,Ka,x,y)+Va[x][y];

//*(Va[x-1][y]-2*Va[x][y]+Va[x+1][y]) + (G)*(Va[x][y-1]-2*Va[x][y]+Va[x][y+1])) +Va[x][y]+IAPP(dt*k,x,y);



            }

        }

        //Border condition

        //  x=0 y=0
        V[0][0]=V[1][0];
        // x=nx-1 y=0
        V[nx-1][0]=V[nx-2][0];
        //  x=0 y=ny-1
        V[0][ny-1]=V[0][ny-2];
        // x=nx-1 y=ny-1
        V[nx-1][ny-1]= V[nx-2][ny-2];

        //y=0 e y=ny para nx-1>x>0
        for(int x=0; x<nx-1; x++)
        {

            //y=0
            V[x][0]=V[x][1];
            V[x][ny-1]= V[x][ny-2];

        }

        //x=0 e x=nx para ny-1>y>y
        for(int y=0; y<=ny-1; y++)
        {

            //x=0
            V[0][y]=V[1][y];
            //x=nx-1
            V[nx-1][y]=V[nx-2][y];
        }



        for(int x=0; x<=nx-1; x++)
        {
            for(int y=0; y<=nx-1; y++)
            {
                K[x][y]=g(V[x][y],Ka[x][y]);

            }

        }


        //copia do vetor

        for (int x=0; x<=nx-1; x++)
        {
            for (int y=0; y<=ny-1; y++)
            {

                if(k%1==0)
                    fprintf(arquivo,"%.3f %.1f %.1f %f \n",dt*k,x*dx,y*dy,V[x][y]);

                if(x==57 && y==57 )

                    fprintf(arquivo2,"%f %f %f %f \n",dt*k,V[x][y],difusao(Va,Ka,x,y),dt*f(Va[x][y],Ka[x][y]));

                Va[x][y] = V[x][y];
                Ka[x][y] =K[x][y];
            }
            fprintf(arquivo," \n");

        }


        if(k%1==0)
            fprintf(arquivo," \n\n");

    }


    fclose(arquivo);





    system("gnuplot -p \"script.txt\" ");
    return 0;
}
