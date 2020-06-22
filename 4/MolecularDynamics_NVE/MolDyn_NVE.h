/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
using namespace std;

//parameters, observables
const int n_props=5;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp,stima_pres;
string Phase;
double walker[n_props];
// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//restart
bool restart;
void restartMode(void);
void standardMode(void);
void ConfRestart(void);

//rescaling
bool rescaling;
double dispv2, actualT,s_f;
double x0old[m_part],y0old[m_part],z0old[m_part];
void rescalingVelocities(double targetT);
void standardVelocities();
double targetT;

//average properties
double blk_av[n_props],blk_norm,accepted,attempted;
double glob_av[n_props],glob_av2[n_props];
double stima_potf,stima_presf,err_epot,err_pres,stima_etotf,err_etot,stima_tempf,err_temp,stima_ekinf,err_ekin;
int b;
int imes=10;
//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
