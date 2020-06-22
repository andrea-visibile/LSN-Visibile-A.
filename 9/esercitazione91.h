

#ifndef esercitazione91
#define esercitazione91

using namespace std;


//input
string config;
const int n=32;
double pm, pc, depth,l1;
const int pop=800;

//functions
void Input();
void BuildingConfig();

//
int seed[4];
double best_m;
double x[pop][n],y[pop][n],e[pop];
double x_new[pop][n],y_new[pop][n];
Random rnd;
void Mutation(int);
void BuildingConfig();
void PairPermutation(int);
void Generation();
void clean();
void MultiplePermutation(int);
void Crossover();
void ContinuosPermutation(int);
void Inversion(int);
double Evaluation(int);
void Print(string c);
void bestl2();
double minl2=9999;

#endif



