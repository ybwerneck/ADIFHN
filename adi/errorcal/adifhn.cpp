#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <bits/stdc++.h>

double dx,dy;
double G;
double X;
double dt;
int Segundos;
int k;
double Y;
double L;
int tam;
double Cm;
double y,a,e;
bool fliped;
int n;




void dump_vectore(double *u,int t)
{
    int j;

    printf("[");
    for(j=0; j<tam; j++)


        printf("(i %d j% d - %f \n ",t,j,u[j]);

    printf("]");


}
void dump_vector(double *u,int tam)
{
    int j;
    printf("[");

    for(j=0; j<tam; j++)
        printf(" %f ",u[j]);


    printf("]");
}
void dump_matrix(double **u,int tamx,int tamy)
{
    int i,j;
    for(i=0; i<tamx; i++)
    {
        printf("[");
        for(j=0; j<tamy; j++)
            printf(" %f ",u[i][j]);
        printf("] \n");
    }
}

void copyMatrix(double **origem,double **destino,int tamx,int tamy)
{
    int i,j;
    for(i=0; i<tam; i++)
    {
        for(j=0; j<tam; j++)
        {
            destino[i][j]=origem[i][j];

        }
    }

}
void flipmatrix(double **u,int tamx,int tamy)
{


    int i,j;
    double** L=(double**) calloc(sizeof(double*), tam);
    for(i=0; i<tam; i++)
    {
        L[i]=(double*) calloc(sizeof(double), tam);
        for(j=0; j<tam; j++)
        {
            L[i][j]=u[j][i];
        }
    }
    copyMatrix(L,u,tam,tam);
    for(i=0; i<tam; i++)
    {

    free(L[i]);
    }
    free(L);
}

void dump_matrix(double **u)
{
    dump_matrix(u,tam,tam);
}
double** initU()
{
    int i,j;
    double** U=(double**) calloc(sizeof(double*), tam);
    for(i=0; i<tam; i++)
    {
        U[i]=(double*) calloc(sizeof(double), tam);
        for(j=0; j<tam; j++)
        {
            U[i][j]=0.0;

            int r=tam/10,c=(tam/2)-1;
            if(((i-c)*(i-c) + (j-c)*(j-c)) < r*r )
            {
                U[i][j]=1;

            }

        }
    }
    return U;
}

double** initG()
{

    int i,j;
    double** G=(double**) calloc(sizeof(double*), tam);
    for(i=0; i<tam; i++)
    {
        G[i]=(double*) calloc(sizeof(double), tam);
        for(j=0; j<tam; j++)
        {
            G[i][j]=0;
        }
    }
    return G;
}
double** iniciarMatrizRHS()
{
    double** matrix =(double**) malloc(sizeof(double*) * tam);
    int i,j;
    matrix[0]=(double*) calloc(sizeof(double), tam);
    matrix[1]=(double*) calloc(sizeof(double), tam);
    matrix[2]=(double*) calloc(sizeof(double), tam);

    for(i=0; i<tam; i++)
    {


        matrix[0][i]=1+(2*Y);
        matrix[1][i]= -(Y);
        matrix[2][i]= -(Y);
    }

    matrix[0][0]=1+Y;
    matrix[0][tam-1]=1+Y;
    matrix[1][0]=0.0;
    matrix[1][tam-1]=-Y;
    matrix[2][0]=-Y;
    matrix[2][tam-1]=0.0;


    return matrix;
}
double aproxU(double U)
{
    return U;
}
double Ic(double U,double g)
{
    return (U*(1 - U)*(U - a) - g)/e;
}
double R(double U,double g)
{
    return Ic(aproxU(U),g) ;

}
float Iapp(double x, double y, double t,double V)
{
    if((t>=0.1 && t<=0.3 ) && (( x>2.8 && x<3.2 )) && (( y>2.8 && y<3.2 )))
    {
        printf("%f %f \n",x,y);
        return  1;
    }

    else
        return V;

}
double calcularFX(double **U,double g,int x,int y )

{

    if (x!=0 && x!=tam-1)
        return U[x][y]*(1 - 2*Y) + Y*U[x+1][y] + Y*U[x-1][y] + R(aproxU(U[x][y]),g)*dt;
    else if(x==0)
        return U[x][y]*(1 - Y) +  Y*U[x+1][y] + R(aproxU(U[x][y]),g)*dt;

    else return U[x][y]*(1 - Y) + Y*U[x-1][y] + R(aproxU(U[x][y]),g)*dt;

}
void getVetorFx(double **U,double **g,int x)
{
    int i;
    for(i=0; i<tam; i++)
       U[x][i]= calcularFX(U,g[x][i],x,i);

    }

void updateG(double **V,double **g,int tam)
{
    int i,j;
    for(i=0; i<tam; i++)
        for(j=0; j<tam; j++)
            g[i][j]= g[i][j] + dt *0.5*(V[i][j]-(g[i][j]*y));

}
void solve_tridiagonal(double x[], const size_t N, const double a[], const  double b[], const  double c[])
{
    size_t n;

    double * const cprime = (double *) malloc(sizeof(double) * N);
    cprime[0] = c[0] / b[0];
    x[0] = x[0] / b[0];

    for (n = 1; n < N; n++)
    {
        double m = 1.0f / (b[n] - a[n] * cprime[n - 1]);
        cprime[n] = c[n] * m;
        x[n] = (x[n] - a[n] * x[n - 1]) * m;
    }

    for (n = N - 1; -n-- > 0; )
        x[n] = x[n] - cprime[n] * x[n + 1];

    free(cprime);
}
char * adifhn(int z, int kt, char c[])
{

    FILE *arquivo;
    FILE *arquivo2;

    char filename[80];
    char filenameb[80];

    if(kt!=-1)
    {

        sprintf(filename,"result_z%d.txt",z);
        sprintf(filenameb,"result_sc_z%d.txt",z);

    }
    else
    {
        sprintf(filename,"result.txt");
        sprintf(filenameb,"result_sc.txt");

    }
    arquivo= fopen(filename,"w");
    arquivo2= fopen(filenameb,"w");
    if(kt!=-1)
        fprintf(arquivo,"%f %f",dt,dx);


    dx=0.1/z,dy=0.1/z;
    L=3;

    tam=L/dx + 2;

    G=25;
    X=100;
    dt=0.01/(z*z);
    Segundos=1;
    k=Segundos/dt;
    Y = (G/X)*((dt)/(2*dx*dx));
    Cm=0.0;
    y=0.5, a=0.1,e=0.01;
    double** u;
    double** g;
    double** rhs=iniciarMatrizRHS();
    u=initU();
    g=initG();

    fliped=false;
    int n;


    for(n=0; n<k; n++)
    {
        system("cls");
        printf("DT: %f DX:%f \n Quadro %d \ %d   completed %f",dt,dx,n,k,((double)n/k)*100);

        for(int i=0; i<tam; i ++)
        {
            getVetorFx(u,g,i);
            solve_tridiagonal(u[i], tam, rhs[1],rhs[0],rhs[2]);
        }

        updateG(u,g,tam);

        fliped=true;
        flipmatrix(u,tam,tam);
        flipmatrix(g,tam,tam);
        for(int i=0; i<tam; i ++)
        {
            getVetorFx(u,g,i);

            solve_tridiagonal(u[i], tam, rhs[1],rhs[0],rhs[2]);
        }
        flipmatrix(u,tam,tam);
        flipmatrix(g,tam,tam);

        fliped=false;

        updateG(u,g,tam);

        for(int v=0; v<tam; v++)
        {
            for(int j=0; j<tam; j++)
            {

                double x=v*dx;
                double y=j*dx;
                double t=n*dt;
               if(kt==-1 || kt ==n)
                {
                    if(kt==-1)
                        fprintf(arquivo,"%.6f %.6f %.6f %f \n",t,x,y,u[v][j]);
                    else
                        fprintf(arquivo,"%f ",u[v][j]);
                            }
            }
            if(kt==-1)
                fprintf(arquivo,"\n");
        }

        if(kt==-1)
        {
           fprintf(arquivo,"\n\n");

        }
        if (kt==n)
            break;

    }

    free(u);
    free(g);

    fclose(arquivo);
    fclose(arquivo2);

    if(kt==-1)
    {

        char gnuplotparam[80];
        sprintf(gnuplotparam,"gnuplot  -e dt=%f -e dx=%f -p  \"script.txt\" %f",dt,dx);

        system(gnuplotparam);

        remove(filename);
        remove(filenameb);
    }
    else
    {
        int l;
        for(l=0; (l<80 || filename[l]=='/0'); l++)
            c[l]=filename[l];
        printf("%s foi salvo!",filename);
    }
    //script 1 anima��o da superficie
    //script 2 grafico do ponto
    //script 3 gif da superficies


}
