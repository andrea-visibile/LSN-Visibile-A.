#ifndef esercitazione51
#define esercitazione51
using namespace std;


//classe random
Random rnd;


//variabili
int acc=0;
double delta;
double x[3];
double x_0,y_0,z_0;
string orbital;
string function;
//quantità per medie a blocchi
const int M=1000000;
const int N=M/2000;
const int L=2000;
double av1[N],av2[N],cum1[N],cum2[N],err[N];
double pos[N][3];
//funzioni

//densità di prob
double hidrogen100(double,double);
double hidrogen210(double,double);

//puntatori
double (*f)(double,double);
double (*t)(double);

//metropolis
double metropolis(double[]);
double uniform(double);
double gaussian(double);
void settings();

//media a blocchi
double error(double[],double[],int);
#endif
