#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include <string>
#include "esercitazione81.h"
#include <iomanip>
using namespace std;

//Esercizio 8.1
//


int main(){
	  rnd.LoadSeed();
	  settings(); //Inizialization
	for(int i=0; i<par;i++){
		for(int j=0; j<par;j++){
		sigma=sigmapar[i];
		mu=mupar[j];
		cout <<"---------------------------------" << endl;
		cout << "Mu: " << mu << "  " << "sigma: " << sigma << endl;
		x[0]=x_0;
			for(int iblk=1; iblk <= nblock; ++iblk){
				   Reset(iblk);   //Reset block averages
				   for(int istep=1; istep <= L; ++istep)
					    {
					      metropolis(x[0]);
					      Measure(x[0]);
					      Accumulate(); //Update block averages
					    }
				   Averages(iblk);   //Print results for current block
			
		  	}
		//cout << glob_av[0]/(double)nblock << endl;
		best(nblock);
		}
		
		
	}
	//Salvo il migliore
	mu=bestmu;
	sigma=bestsigma;
	x[0]=x_0;
	for(int iblk=1; iblk <= nblock; ++iblk){
				   Reset(iblk);   //Reset block averages
				   for(int istep=1; istep <= L; ++istep)
					    {
					      metropolis(x[0]);
					      Measure(x[0]);
					      Accumulate(); //Update block averages
					      Print_x();
						//cout << H << endl;	
					    }
				   Averages(iblk);   //Print results for current block
				
				   Print(iblk);
				   
		  	}

	
}

/////////////////////////////////*//////////////////////////////////////////////
//funzioni 

void best(int nblock){
	ofstream out;
	
	//cout <<glob_av[0]/(double)nblock<< " " << amin << endl;
	if ( glob_av[0]/(double)nblock < amin) {
		bestmu=mu;
		bestsigma=sigma;
		out.open("bestpar.dat");
		amin=glob_av[0]/(double)nblock;
		//cout << "printing on file" << endl;
		out << mu << setw(12) << sigma << endl;
		out.close();
		}

}


double Accumulate(){

	blk_av[0] = blk_av[0] + H;
	blk_norm = blk_norm + 1.0;
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
 	glob_av[0] = 0;
	glob_av2[0] = 0;
   }


   blk_av[0] = 0;
   blk_norm = 0;
   attempted = 0;
   acc = 0;
}

double Measure(double x){
	H=-0.5 * der(x,mu,sigma)/doublegauss(x,mu,sigma)+V(x);
	return H;
}

double der(double x, double mu, double sigma){
	double m1,m2,exp1,exp2;
	m1 = pow(((x-mu)/sigma),2);
	m2 = pow(((x+mu)/sigma),2);
	exp1 = exp(-m1/2.);
	exp2 = exp(-m2/2.);
	return (exp1*(m1-1.) + exp2*(m2-1.))/(sigma * sigma);

}

//calcolare l'errore nella media a blocchi
double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

//distribuzione uniforme
double uniform(double x){

	return (x  + (rnd.Rannyu(-delta,delta)));

}

//distribuzione gaussiana
double gaussian(double x){

	return rnd.Gauss(x,delta);

}

//funzione di prob
double doublegauss(double x,double mu, double sigma){
	double d,m1,m2;
	m1=pow((x-mu),2);
	m2=pow((x+mu),2);
	d = 2*sigma*sigma;
	return exp(-m1/d)+exp(-m2/d);
}

//modulo quadro della funzione di prob
double doublegauss2(double x, double mu, double sigma){
	return pow(doublegauss(x,mu,sigma),2);
}

double V(double x){
	return pow(x,4)-(5./2.)*x*x;
}


//Metropolis!
void metropolis(double xo){
        double y;
	double p,r,px,py;
	
	y= t(xo);

	

	px=f(xo,mu,sigma);
	py=f(y,mu,sigma);
	p=min(1.,py/px);
	

	r=rnd.Rannyu();
	
	if (r<=p){
		acc++;
		x[0]=y;
	}
	attempted++;

}
void Averages(int iblk){
	
	//cout << "Block number " << iblk << endl;
	//cout << "Acceptance rate " << acc*1./attempted << endl << endl;


	stima_H = blk_av[0]/(1.*blk_norm);
	glob_av[0] += stima_H;
	glob_av2[0] += stima_H*stima_H;
	err_H= Error(glob_av[0],glob_av2[0],iblk);
}


void Print(int iblk){
	ofstream out;
	out.open("H.dat",ios::app);
	out << setw(12) << iblk <<  setw(12) << stima_H << setw(12) << glob_av[0]/(double)iblk << setw(12) << err_H << endl;
}

void Print_x(){
	ofstream out;
	out.open("x.dat",ios::app);
	out << x[0] << endl;
}
//Funzione per evitare di dover cambiare ogni volta tutti i parametri
//settings.C contiene informazioni sull'orbitale e sulla funzione da usare in metropolis
void settings(){
	ifstream in;
	in.open("settings.C");
	cout << "Metropolis algorithm to solve a 1D quantum problem" << endl;
	


	cout << "Reading input from settings.C" << endl;
	in >> function;
	cout << "You are using a uniform probability" << endl;
	in >> delta;
	cout << "The step of the simulation is " << delta << endl;
	in >> x_0;
	cout << "The starting point is " << x_0 << endl;
	in >> nsteps;
	cout << "Simulation with " << nsteps << " step" << endl;
	in >> nblock;
	cout << "Simulation divided in  " << nblock << " blocks" << endl;
	L=nsteps/nblock;

	


	if (function == "uniform") t=&uniform;
	if (function =="gaussian") t=&gaussian;

	f=&doublegauss2;

	//grid of parameters
	for (int i=0; i<par;i++){
		sigmapar[i] = 0.5+i*0.002;
		mupar[i] = 0.6+i*0.002;
	}

}

