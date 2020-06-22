



using namespace std;


//input
string config;
int size,seed_mpi;
const int n=32;
double pm, pc, depth,l1;
const int pop=800;
int argmin;
//functions
void Input();
void BuildingConfig();
void Shuffle_config_mpi();
void Transmission();
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
void Print();
void bestl2();
double minl2=9999;

class city{
	double x;
	double y;
	};

class path{
	city p[32];
	double eval;
	};




