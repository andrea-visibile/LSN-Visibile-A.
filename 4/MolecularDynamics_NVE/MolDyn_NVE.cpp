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

   for(int iblk=1; iblk <= nstep/b; ++iblk) //Simulation
  {
    if(iblk%iprint == 0) cout << "Number of blocks: " << iblk << endl;
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= b; ++istep)
    {
	
        Move();
	if(istep%20 == 0){
	     Measure();     
	     Accumulate();
	     nconf += 1;
	}
    }
    Averages(iblk);   //Print results for current block
  }
  ConfRestart();
  ConfFinal();         //Write final configuration to restart

  return 0;

}
////////////////////////////////////////*********************************************///////////////////////////////////////////////////

//funzioni
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


void Accumulate(void){

walker[0]=stima_temp;
walker[1]=stima_kin;
walker[2]=stima_pot;
walker[3]=stima_etot;
walker[4]=stima_pres;

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream  Epot, Pres;
   const int wd=12;
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
		  Pres_out.open("Argon/"+Phase+"/ave_pres.out",ios::app);
		  Pres_err_out.open("Argon/"+Phase+"/ave_err_pres.out",ios::app);
	}
    //cout << "Block number " << iblk << endl;
	
    stima_tempf = blk_av[0]/blk_norm ; //temp
    glob_av[0] += stima_tempf;
    glob_av2[0] += stima_tempf*stima_tempf;
    err_temp=Error(glob_av[0],glob_av2[0],iblk);


    Temp_out << glob_av[0]/(double)iblk << endl;
    Temp_err_out << err_temp << endl;
	
    stima_ekinf = blk_av[1]/blk_norm ; //kinetic energy
    glob_av[1] += stima_ekinf;
    glob_av2[1] += stima_ekinf*stima_ekinf;
    err_ekin=Error(glob_av[1],glob_av2[1],iblk);

    Ekin_out << glob_av[1]/(double)iblk << endl;
    Ekin_err_out << err_ekin << endl;
    
    stima_potf = blk_av[2]/blk_norm ; //Potential energy
    glob_av[2] += stima_potf;
    glob_av2[2] += stima_potf*stima_potf;
    err_epot=Error(glob_av[2],glob_av2[2],iblk);

    Epot_out << glob_av[2]/(double)iblk << endl;
    Epot_err_out << err_epot << endl;
	
    stima_etotf = blk_av[3]/blk_norm ; //total energy
    glob_av[3] += stima_etotf;
    glob_av2[3] += stima_etotf*stima_etotf;
    err_etot=Error(glob_av[3],glob_av2[3],iblk);

    Etot_out << glob_av[3]/(double)iblk << endl;
    Etot_err_out << err_etot << endl;
        

    stima_presf = (blk_av[4]/blk_norm); //Pressure
    glob_av[4] += stima_presf;
    glob_av2[4] += stima_presf*stima_presf;
    err_pres=Error(glob_av[4],glob_av2[4],iblk);

    Pres_out << glob_av[4]/(double)iblk << endl;
    Pres_err_out << err_pres << endl;

}


double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
    stima_pres = stima_temp*rho + w/(vol*3.);
	
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
