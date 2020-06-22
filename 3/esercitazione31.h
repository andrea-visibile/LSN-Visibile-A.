#ifndef esercitazione31
#define esercitazione31
using namespace std;

Random rnd;

const int M = 10000; //steps
const int N=100; //blocks
int L= M/N;


//parametri
double S_0 = 100; //prezzo al tempo iniziale
double T = 1;
double K = 100; //strike price
double r = 0.1; //interest rate 
double vol = 0.25; //volatility
double av1[N],av4[N],av3[N];
double av2[N];
double cum1[N], cum3[N];
double cum2[N], cum4[N];
double err[N], err2[N];


//punto a (direct sampling)
double S;
double throws[M];
double throws2[M];
double error(double av1[], double av2[], int n);
#endif
