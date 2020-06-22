
#ifndef esercitazione22
#define esercitazione22

using namespace std;

Random rnd;

//numero di tiri totali
const int n = 10000;
const int step = 100;
const int N=100;
//assi
double x[3];

//altro
double r[n];
int ax;


//punto b
double theta, phi;
double x2[3];

//vettori in cui salvare le medie sui blocchi
double av1[step][n],av4[step][n],av3[step][n],av2[step][n];
double sum = 0;
double sum2 = 0;
double cum1[N], cum3[N];
double cum2[N], cum4[N];
double err[N], err2[N];
double error(double av1[], double av2[], int n);

#endif
