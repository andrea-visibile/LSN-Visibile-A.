/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
#include <string>

using namespace std;

int main()


{ 
  int soglia=30; //Salta le prime 30 misure
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      if(iblk>soglia){
		//iblk=iblk-soglia;
	      //cout << "measuring block: " << iblk << endl;
	      Measure();
	      Accumulate(); //Update block averages
		//iblk=iblk+soglia;
	}
    }
   if(iblk>soglia){
    Averages(iblk-soglia);}   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput >> restart;
  if (restart==1) cout <<"Reading old config from a previous run" << endl;
  ReadInput >> prova;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
if (restart==0){


  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
 
  }
}
else{
  ifstream out;
  out.open("config.final");
  if (!out.is_open()) cerr << "PROBLEM: Unable to open config.final." <<endl;
  cout << "reading spins" << endl;
  for (int i=0; i<nspin; ++i){
  	out >> s[i];
	
	}
  out.close();


}


//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle 
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
	
 	metropolis(o);
    }
    else //Gibbs sampling
    {
	
	gibbs(o);
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
 
  double u = 0.0, si = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     si += s[i];

  }
  // walker le grandezze necessarie a calcolare le proprietà termodinamiche (gli indici sono 0,1,2,3)
  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[im] = si;
  walker[ix] = pow(si,2);
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Cal, Magn, Susx;
   const int wd=18;
   string oo;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    if (metro==1) oo="_metro"; //per salvare in file diversi metropolis e gibbs
	else oo="_gibbs";
    if(prova=="prova"){//mi serve per i primi grafici sull'equilibrazione e per fare prove senza salvare nella cartella temp
	    Ene.open("a.0",ios::app);

	    Susx.open("b.0",ios::app);
	    Cal.open("c.0",ios::app);
	    Magn.open("d.0",ios::app);
	}
    else{
		if(h==0){
		    Ene.open("temp/output"+to_string(temp)+oo+".ene.0",ios::app); //per evitare di sovrascrivere magnetizzazione quando faccio a h=0
		    Susx.open("temp/output"+to_string(temp)+oo+".Susx.0",ios::app);
		    Cal.open("temp/output"+to_string(temp)+oo+".Cal.0",ios::app);
			}
		else{
	    	Magn.open("temp/output"+to_string(temp)+oo+".Magn.0",ios::app);
		}
	}
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();


    stima_c = beta*beta*((blk_av[ic]/blk_norm)-(blk_av[iu]/blk_norm)*(blk_av[iu]/blk_norm)); 
    stima_c/=(double)nspin;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Cal << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Cal.close();

 
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Mag
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Magn << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Magn.close();

    
    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Susx << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Susx.close();



    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
	
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void metropolis(int o){
	
	double p,dE,r;
	dE= -2.*Boltzmann(s[o],o); //normale metropolis

	p=min(1.,exp(-dE*beta));
	

	r=rnd.Rannyu();
	if (r<=p){
		accepted++;
		s[o] = -s[o];
	}

	attempted++;
}

void gibbs(int o){
      //Algoritmo di Gibbs (heat-Bath)
      double dE,q,p,r;

      dE=-2.*Boltzmann(1,o); //uso il meno
      q=exp(-beta*dE);
      p=1./(1.+q); 

      r=rnd.Rannyu();

      if(r<=p)  s[o]=1.;//quindi più
      else s[o]=-1.; //anche in caso else muove lo spin: l'accettanza è per forza 1
      attempted ++; 
      accepted ++;




}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
