
#ifndef esercitazione81
#define esercitazione81
using namespace std;


//classe random
Random rnd;


//variabili
int acc,attempted;
double delta;
double x[3];
double x_0,y_0,z_0;
string orbital;
string function;
int nblock, nsteps, L;
double amin=9999;
double bestmu, bestsigma;
double mu=1;
double sigma=1;
double H, stima_H,err_H;
const int m_props=1000;
const int par=150;
//funzioni

//densità di prob
double doublegauss(double, double, double);
double doublegauss2(double, double, double);
double V(double);
double der(double x, double mu, double sigma);

//blocchi
double Measure(double x);
double Accumulate();
void Reset(int);
void Print(int);
void Print_x();
void Averages(int);
double blk_norm;
double blk_av[m_props], glob_av[m_props], glob_av2[m_props];
//puntatori
double (*f)(double,double,double);
double (*t)(double);

//metropolis
void metropolis(double);
double uniform(double);
double gaussian(double);
void settings();
//parametri
double mupar[par];
double sigmapar[par];
void best(int);
//media a blocchi
double error(double[],double[],int);

#endif
