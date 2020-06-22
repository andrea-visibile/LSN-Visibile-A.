
#ifndef esercitazione31
#define esercitazione31

using namespace std;

// variabili necessarie per l'esperimento di Buffon
double l = 4; //lunghezza sbarra
double d = 5; //distanza fra le linee
double a,y,b; 
int sum = 0;

//variabili necessarie per la media a blocchi
int M=10000;
const int N=100;
int L= M/N;
double av1[N],av2[N];
double cum1[N], cum2[N],err[N];

//funzioni
double pi(int N_hit, int N_throws,double D, double L);
int buffon(double d, double L, double a, double y);
double error(double av1[], double av2[], int n);

#endif
