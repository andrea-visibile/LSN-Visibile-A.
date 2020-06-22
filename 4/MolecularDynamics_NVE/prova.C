/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <string>

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1;
  int n = 0;
  double sum[m_props];

  for (int i=0; i<m_props; i++){
		sum[i]=0;	
		}
  //Vettori	
  double av_etot[nstep/b], av2_etot[nstep/b]; 
  double av_temp[nstep/b], av2_temp[nstep/b];
  double av_epot[nstep/b], av2_epot[nstep/b];
  double av_ekin[nstep/b], av2_ekin[nstep/b];
  double av_pres[nstep/b], av2_pres[nstep/b];


  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%imes == 0){ //ogni tot misura
        Measure();     //Properties measurement 
        nconf += 1;
	sum[0] += stima_pot;
	sum[1] += stima_kin;
	sum[2] += stima_temp;
	sum[3] += stima_etot;
	sum[4] += stima_pres;
     }
     if(istep%b == 0){ //b è il numero di step per blocco, voglio che resetti tutto e calcoli la media del blocco
	av_etot[n]=sum[3]/(b/imes); 
	av2_etot[n]=pow(av_etot[n],2);
	av_epot[n]=sum[0]/(b/imes); 
	av2_epot[n]=pow(av_epot[n],2);
	av_ekin[n]=sum[1]/(b/imes); 
	av2_ekin[n]=pow(av_ekin[n],2);
	av_temp[n]=sum[2]/(b/imes); 
	av2_temp[n]=pow(av_temp[n],2);
	av_pres[n]=sum[4]/(b/imes); 
	av2_pres[n]=pow(av_pres[n],2);
	n++;
	for (int y=0 ; y<4 ; y++){
		sum[y]=0;	
		}
     }
  }
  MeasureAverages(av_etot,  av2_etot, av_temp, av2_temp, av_epot, av2_epot, av_ekin, av2_ekin,av_pres, av2_pres); //valori medi
  ConfRestart();
  ConfFinal();         //Write final configuration to restart

  return 0;


////////////////////////////////////////*********************************************///////////////////////////////////////////////////

//funzioni


}
// funzione per calcolare l'errore con il metodo a blocchi
double error(double av1[], double av2[], int n){
	double err;
	if (n==0){
		return 0;
		}
	else {
		err=av2[n]-pow(av1[n],2);
		return sqrt(err/n);
		}
	}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input
  ReadInput >> Phase;
  cout << "You are studying a " << Phase << " phase of Argon! " << endl << endl << endl;
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;
  ReadInput >> rescaling;
  ReadInput >> targetT;
  ReadInput >> b;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  cout << "This is " << restart << " mode" << endl;
  cout << "This is " << rescaling << " mode" << endl;
  cout << "The rescaling temperature has been set to " << targetT << endl;
  cout << "Number of step in a block (please choose a multiple of 10): " << b << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables



  if(restart == false){
	standardMode();
	}
  else{
	cout <<"entering restard mode" << endl;
	restartMode();
	}

  if(rescaling == false ){
	standardVelocities();
	}
  else{
	rescalingVelocities(targetT);
	}
	
	


}

//Read initial configuration
void standardMode(){
	ifstream ReadConf;
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
		}
	ReadConf.close();
}

//Read supplementary configuration for restart mode 
void restartMode(){
	ifstream ReadConf;
	cout <<"Restart mode! " << endl;	
	cout << "Read supplementary configuration from file old.0 " << endl << endl;
	ReadConf.open("old.0");

	//leggo la configurazione old (ovvero quella prima della xi)
	if(!ReadConf.is_open()) {
	    cerr << "It seems the file is missing!" << endl;
	    exit(1);
	  }
	for (int i=0; i<npart; ++i){
		ReadConf >> xold[i] >> yold[i] >> zold[i];
		xold[i] = xold[i] * box;
		yold[i] = yold[i] * box;
		zold[i] = zold[i] * box;
		}
	ReadConf.close();
	
	//leggo la configurazione ultima a cui è stato portato xi nella precedente run
	ReadConf.open("old.final");

	if(!ReadConf.is_open()) {
	    cerr << "It seems the file is missing! The algorithm required that you create a starting point via standardMode simulation!" << endl;
	    exit(1);
	  }
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
		}
	ReadConf.close();
}

void standardVelocities(){
	//Prepare initial velocities
	  
	   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	   double sumv[3] = {0.0, 0.0, 0.0};
	   for (int i=0; i<npart; ++i){
	     vx[i] = rand()/double(RAND_MAX) - 0.5;
	     vy[i] = rand()/double(RAND_MAX) - 0.5;
	     vz[i] = rand()/double(RAND_MAX) - 0.5;

	     sumv[0] += vx[i];
	     sumv[1] += vy[i];
	     sumv[2] += vz[i];
	   }
	   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	   double sumv2 = 0.0, fs;
	   for (int i=0; i<npart; ++i){
	     vx[i] = vx[i] - sumv[0];
	     vy[i] = vy[i] - sumv[1];
	     vz[i] = vz[i] - sumv[2];

	     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	   }
	   sumv2 /= (double)npart;

	   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	   for (int i=0; i<npart; ++i){
	     vx[i] *= fs;
	     vy[i] *= fs;
	     vz[i] *= fs;
		
	     xold[i] = Pbc(x[i] - vx[i] * delta);
	     yold[i] = Pbc(y[i] - vy[i] * delta);
	     zold[i] = Pbc(z[i] - vz[i] * delta);
	   }
	   
}


void rescalingVelocities(double targetT){
	cout << "Rescaling velocities!" << endl;
	Move(); //faccio un passo con l'algoritmo
	

	double sumv2 = 0.0; //calcolo le velocità
	for (int i=0; i<npart; ++i){

	     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	   }
	sumv2 /= (double)npart;
	s_f = sqrt(3 * targetT / sumv2);
 
	for(int i=0; i<npart; ++i){ //le riscalo
	    vx[i] *= s_f;
	    vy[i] *= s_f;
	    vz[i] *= s_f;
 
  	    double x_i=xold[i]; 
	    double y_i=yold[i];
	    double z_i=zold[i];
	    
	    xold[i] = Pbc(xold[i] - (vx[i])*delta); //da xi torno indietro  e mi calcolo un nuovo old tramite le velocità
	    yold[i] = Pbc(yold[i] - (vy[i])*delta);
	    zold[i] = Pbc(zold[i] - (vz[i])*delta);
	     
	    x[i]=x_i; //assegno xold a x, poichè ho già fatto un passo
	    y[i]=y_i;
	    z[i]=z_i;
  
	    }

	
	

}

void MeasureAverages(double av_etot[], double av2_etot[],double av_temp[], double av2_temp[], double av_epot[], double av2_epot[],double av_ekin[], double av2_ekin[], double av_pres[],double av2_pres[]){


  double cum_pot[b], cum2_pot[b], cum_tot[b], cum2_tot[b], cum_kin[b], cum2_kin[b], cum_temp[b], cum2_temp[b],cum_pres[b],cum2_pres[b]; 
  double err_etot[nstep/b], err_epot[nstep/b], err_ekin[nstep/b], err_temp[nstep/b], err_pres[nstep/b];


  for (int i=0; i<b; i++){ //azzero i vettori con la somma
		cum_pot[i]=0;
		cum2_pot[i]=0;	
		cum_tot[i]=0;
		cum2_tot[i]=0;	
		cum_kin[i]=0;
		cum2_kin[i]=0;
		cum_temp[i]=0;
		cum2_temp[i]=0;
		cum_pres[i]=0;
		cum2_pres[i]=0;
		}
  
  //calcolo la grandezza al variare del numero di blocchi
	for (int i=0; i<(nstep/b); i++){ 
		for (int j=0; j<i+1; j++){ 
			cum_pot[i] += av_epot[j]; 
			cum2_pot[i] += av2_epot[j]; 
			cum_tot[i] += av_etot[j]; 
			cum2_tot[i] += av2_etot[j]; 
			cum_kin[i] += av_ekin[j]; 
			cum2_kin[i] += av2_ekin[j]; 
			cum_temp[i] += av_temp[j]; 
			cum2_temp[i] += av2_temp[j];
			cum_pres[i] += av_pres[j]; 
			cum2_pres[i] += av2_pres[j];  
			}
		cum_pot[i] /= (i+1); //faccio la media di i+1 blocchi
		cum2_pot[i] /= (i+1); 
		err_epot[i]= error( cum_pot, cum2_pot, i); //vettore di errori 

		cum_tot[i] /= (i+1); //faccio la media di i+1 blocchi
		cum2_tot[i] /= (i+1); 
		err_etot[i]= error( cum_tot, cum2_tot, i);

		cum_kin[i] /= (i+1); //faccio la media di i+1 blocchi
		cum2_kin[i] /= (i+1); 
		err_ekin[i]= error( cum_kin, cum2_kin, i);

		cum_temp[i] /= (i+1); //faccio la media di i+1 blocchi
		cum2_temp[i] /= (i+1); 
		err_temp[i]= error( cum_temp, cum2_temp, i);

		cum_pres[i] /= (i+1); //faccio la media di i+1 blocchi
		cum2_pres[i] /= (i+1); 
		err_pres[i]= error( cum_pres, cum2_pres, i);
		}
	//print on file
	ofstream Epot_out, Ekin_out, Etot_out, Temp_out, Pres_out;
  	ofstream Epot_err_out, Ekin_err_out, Etot_err_out, Temp_err_out,Pres_err_out;
	//voglio file diversi per l'argon e per la simulazione semplice necessaria per l'equilibrazione
	//voglio anche che a seconda della fase me le stampi direttamente nelle cartelle
	if(Phase=="0"){
		  Epot_out.open("ave_epot.out",ios::app);
		  Ekin_out.open("ave_ekin.out",ios::app);
		  Etot_out.open("ave_etot.out",ios::app);
		  Temp_out.open("ave_temp.out",ios::app);
		  Epot_err_out.open("ave_err_epot.out",ios::app);
		  Ekin_err_out.open("ave_err_ekin.out",ios::app);
		  Etot_err_out.open("ave_err_etot.out",ios::app);
		  Temp_err_out.open("ave_err_temp.out",ios::app);
		  Pres_out.open("ave_pres.out",ios::app);
		  Pres_err_out.open("ave_err_pres.out",ios::app);
		
		}
	else{
		  Epot_out.open("Argon/"+Phase+"/ave_epot.out",ios::app);
		  Ekin_out.open("Argon/"+Phase+"/ave_ekin.out",ios::app);
		  Etot_out.open("Argon/"+Phase+"/ave_etot.out",ios::app);
		  Temp_out.open("Argon/"+Phase+"/ave_temp.out",ios::app);
		  Epot_err_out.open("Argon/"+Phase+"/ave_err_epot.out",ios::app);
		  Ekin_err_out.open("Argon/"+Phase+"/ave_err_ekin.out",ios::app);
		  Etot_err_out.open("Argon/"+Phase+"/ave_err_etot.out",ios::app);
		  Temp_err_out.open("Argon/"+Phase+"/ave_err_temp.out",ios::app);
		  Pres_out.open("Argon/"+Phase+"ave_pres.out",ios::app);
		  Pres_err_out.open("Argon/"+Phase+"ave_err_pres.out",ios::app);
	}

	for(int i=0; i<(nstep/b); i++){ //stampo su file
		Epot_out << cum_pot[i]<<endl;
		Ekin_out << cum_kin[i]<<endl;
		Etot_out << cum_tot[i]<<endl;
		Temp_out << cum_temp[i]<<endl;
		Epot_err_out << err_epot[i]<<endl;
		Ekin_err_out << err_ekin[i]<<endl;
		Etot_err_out << err_etot[i]<<endl;
		Temp_err_out << err_temp[i]<<endl;
		Pres_err_out << err_pres[i]<<endl;
		Pres_out << cum_pres[i]<<endl;
	}
	Epot_out.close();
	Etot_out.close();
	Ekin_out.close();
	Temp_out.close();
	Epot_err_out.close();
	Etot_err_out.close();
	Ekin_err_out.close();
	Temp_err_out.close();
	Pres_out.close();
	Pres_err_out.close();
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
  
	
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
	 
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij,wij,w;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;
  if(Phase=="0"){
	  Epot.open("output_epot.dat",ios::app);
	  Ekin.open("output_ekin.dat",ios::app);
	  Temp.open("output_temp.dat",ios::app);
	  Etot.open("output_etot.dat",ios::app);
	  Pres.open("output_pres.dat",ios::app);
	}
  else{
 	  Epot.open("Argon/"+Phase+"/output_epot.dat",ios::app);
	  Ekin.open("Argon/"+Phase+"/output_ekin.dat",ios::app);
	  Temp.open("Argon/"+Phase+"/output_temp.dat",ios::app);
	  Etot.open("Argon/"+Phase+"/output_etot.dat",ios::app);
	  Pres.open("Argon/"+Phase+"output_pres.dat",ios::app);
	if(!Epot.is_open()) {
	    cerr << "It seems the file is missing!!!!!!!!!" << endl;
	    exit(1);
	  }
	}


  v = 0.0; //reset observables
  t = 0.0;
  w=0.;
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 48.0/pow(dr,12) - 24.0/pow(dr,6);
//Potential energy
       v += vij;
       w += wij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = stima_temp*rho + w/(vol*3);

    Pres << stima_pres << endl;
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Pres.close();
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
}









void ConfRestart(void){ //Write final configuration
  ofstream Writeold,Writeold2;

  cout << "Print restart configuration to file old.0 e old.final " << endl << endl;
  Writeold.open("old.0");
  Writeold2.open("old.final");
  for (int i=0; i<npart; ++i){
    Writeold2 << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  for (int i=0; i<npart; ++i){
    Writeold << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  Writeold.close();
  Writeold2.close();
  return;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
