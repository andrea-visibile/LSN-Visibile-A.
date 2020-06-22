
#ifndef esercitazione21
#define esercitazione21

using namespace std;

Random rnd;
//numero di tiri totali
const int M=10000;
//numero di blocchi
const int N=100;
//numero di tiri per blocchi
int L= M/N;

//vettori in cui salvare le medie sui blocchi
double av1[N],av4[N],av3[N],av2[N];
double x,y;
double sum = 0;
double sum2 = 0;
double cum1[N], cum3[N];
double cum2[N], cum4[N];
double err[N], err2[N];

//vettori per uniform e importance sampling
double unif[M];
double imp[M]; 
double error(double av1[], double av2[], int n);

#endif
