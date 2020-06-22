

#ifndef esercitazione11
#define esercitazione11

using namespace std;

//random
Random rnd;

//datablock
const int M=1000000;   //numero di tiri totali
const int N=100; //numero di blocchi
int L= M/N; //numero di tiri per blocchi


double error(double av1[], double av2[], int n);
double throws[M];


//chi
int n=10000; //numero di estratti
const int J=100; //numero di chi
int counter[J];
double chi[J];


//vettori in cui salvare le medie sui blocchi
double av1[N],av4[N],av3[N];
double av2[N];
double cum1[N], cumvar1[N];
double cum2[N], cumvar2[N];
double err[N], errvar[N];

#endif
