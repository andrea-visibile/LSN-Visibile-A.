

#ifndef esercitazione101
#define esercitazione101

using namespace std;


//input
string config;
const int n=32;
double temp, dtemp, nstep,l1;

//functions
void Input();
double t,temp0;
void BuildingConfig();

//
int accepted,attempted;
double pm=1;
int seed[4];
double x[n],y[n];
double x_new[n],y_new[n];
Random rnd;
double beta;
void Mutation();
void Metropolis();
void BuildingConfig();
void PairPermutation();
void clean();
void MultiplePermutation();
void ContinuosPermutation();
void Inversion();
double Evaluation(double x[],double y[]);
void AnnelingSchedule();
void Print(string c);
void bestl2();
double minl2=9999;

#endif




