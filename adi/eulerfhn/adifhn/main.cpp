#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>

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
bool fliped=false;
int n;




void dump_vectore(double *u,int t)
{
    int j;

    printf("[");
    for(j=0; j<tam; j++)


        printf("(i %d j% d - %f \n ",t,j,u[j]);

    printf("]\n");


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

            int r=5*1,c=(tam-1)/2;
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

    if (y!=0 && y!=tam-1)
        return U[x][y]*(1 - 2*Y) + Y*U[x][y-1] + Y*U[x][y+1] + R(aproxU(U[x][y]),g)*dt;
    else if(y==0)
        return U[x][y]*(1 - Y) +  Y*U[x][y+1] + R(aproxU(U[x][y]),g)*dt;

    else if(y==tam-1)
        return U[x][y]*(1 - Y) + Y*U[x][y-1] + R(aproxU(U[x][y]),g)*dt;

}
double* getVetorFx(double **U,double **g,int y)
{
    int i;
    double* Fx = ((double*) calloc(sizeof(double), tam));
    for(i=0; i<tam; i++)
        Fx[i]= calcularFX(U,g[i][y],i,y);

    return Fx;
}
double updateG(double **V,double **g,int tam)
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
void adifhn(int o, int z,int l)
{

    FILE *arquivo;
    FILE *arquivo2;
    char filename[80];
    char filenameb[80];
                sprintf(filename,"result.txt",z);
        sprintf(filenameb,"result_sc.txt",z);
    arquivo= fopen(filename,"w");
    arquivo2= fopen(filenameb,"w");
    dx=0.1/z,dy=0.1/z;
        L=3;

    tam=L/dx + 2;

    G=25;
    X=100;
    dt=0.001/(z*z);
    Segundos=2;
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

        for(int i=0; i<tam; i ++)
        {
            u[i] = getVetorFx(u,g,i);
            solve_tridiagonal(u[i], tam, rhs[1],rhs[0],rhs[2]);
        }

        updateG(u,g,tam);

        fliped=true;
        flipmatrix(u,tam,tam);
        flipmatrix(g,tam,tam);

        for(int i=0; i<tam; i ++)
        {
            u[i] = getVetorFx(u,g,i);
            if(i==0 || i==30)
            dump_vectore(u[i],i);
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
                // if(n%10==0)
                fprintf(arquivo,"%.3f %.1f %.1f %f \n",t,x,y,u[v][j]);
                if(j==57 && v==57 )
                    fprintf(arquivo2,"%f %f %f \n",dt*n,g[v][j],u[v][j]);
            }
            //if(n%10==0)
            fprintf(arquivo,"\n");
        }

//            if(n%10==0)
        fprintf(arquivo,"\n\n");
    }

    free(u);
    free(g);
    fclose(arquivo);
    fclose(arquivo2);

    //script 1 animação da superficie
    //script 2 grafico do ponto
    //script 3 gif da superficies

    switch(l)
    {
    case 1:
        system("gnuplot -p  \"script.txt\" ");
        break;
    case 2:
        system("gnuplot -p \"script2.txt\" ");
        break;
    case 3:
        system("gnuplot -p \"script3.txt\" ");
        break;
    }


}

int main()
{

    adifhn(1,1,1);
}
