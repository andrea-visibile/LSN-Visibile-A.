#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "random.h"
#include "esercitazione101.h"
#include <fstream>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
using namespace std;

int main(){
	Input();
	clean();
	BuildingConfig();
	
	for(int i=0; i<(200); i++){
		cout <<"         *******************            " << endl;
		cout <<"              Actual T= "<< temp << endl; 
		//cout <<" L^(2)_old= " << Evaluation(x,y) << endl;
		for(int i=0; i<nstep; i++){
			//cout <<" L^(2)_old= " << Evaluation(x,y) << endl;
			Mutation();	
			Metropolis();
			//cout <<" L^(2)_new= " << Evaluation(x,y) << endl;
		}
		//cout <<" L^(2)_new= " << Evaluation(x,y) << endl;
		bestl2();
		
		AnnelingSchedule();

	}
	Print("config"+config+".final");
	

return 0;
}

///////////////////////////////////**********************************///////////////////////////////////////////

//Funzioni

void Input(){
	ifstream read;

	//Read seed for random numbers
	rnd.LoadSeed();



	read.open("input.dat");
	cout << endl;
	cout << "***********************************************************************" << endl;
	cout << "                  Traveling Salesman Problem                       " << endl<<endl;
	cout << "                      Simulated Anneling                    "<<endl;
	read >> config;
	cout << "***********************************************************************" << endl;
	cout << endl;
	cout << "The cities configuration is a " << config << endl;

	//read >> n;
	read >> nstep;
	read >> temp0;
	read >> dtemp;
	cout << "There are "<< n << " cities" << endl;
	cout << "Initial temperature: " << temp << endl;
	cout << "Decresing temperature every step: " << dtemp<<endl;
	
	beta=1./temp;
	


}


void AnnelingSchedule(){ //come cambio la temperatura e gli step
	temp=temp0*exp(-dtemp *t);
	nstep=nstep +100;
	t=t+0.1;
}



void Metropolis(){
  
  
  //Probability to accept:
  double deltaL=Evaluation(x_new,y_new)-Evaluation(x,y); //metropolis con energia fittizia
  double q=exp(-beta*deltaL);
  double A=min(1.,q);

  double r=rnd.Rannyu();
  if (r<=A){
    	for(int i=0; i<n; i++){
		x[i]=x_new[i];
		y[i]=y_new[i];
	
	}
    accepted++;
  }
  attempted++;


}



void bestl2(){ //stampo passaggio per passaggio

		
	ofstream out;

	out.open("best"+config+".dat",ios::app);
	out << Evaluation(x,y) << endl;
	out.close();




}


void clean(){ //pulisco il file
	ofstream out;
	out.open("best"+config+".dat");
		out <<"";
	out.close();
	


}

void BuildingConfig(){ //costruisco le città

	if(config == "circle") {
		for(int i=0; i<n; i++){
			x[i] = rnd.Rannyu(-1,1);
			y[i] = pow((-1),i)*sqrt(1-x[i]*x[i]);
			
 			}
		}
	if(config == "square") {
		for(int i=0; i<n; i++){
			x[i] = rnd.Rannyu(0,1);
			y[i] = rnd.Rannyu(0,1);
			}  
		}


}

double Evaluation(double x[],double y[]){ // l2
	l1=0;
	for(int i=0; i<n-1;i++){
		l1 += abs(sqrt(pow((x[i]-x[i+1]),2)+pow((y[i]-y[i+1]),2)));
		//cout << l1 << endl;
		}
	l1 += abs(sqrt(pow((x[n-1]-x[0]),2)+pow((y[n-1]-y[0]),2)));	
	return l1;
}












void Mutation(){ // da mantenere fermo il primo elemento ancora

	double a,r;
	a=rnd.Rannyu(0,1);
	if (a < pm){ r=int(rnd.Rannyu()*4);}
	
	if(r == 0){PairPermutation();}
	if(r == 1){MultiplePermutation();}
	if(r == 2){ContinuosPermutation();}
	if(r == 3){Inversion();}
}

void PairPermutation(){ //primo elemento non permuta

	int r1,r2;
	double appo;
	r1 = int(rnd.Rannyu(1,n-1));
	r2 = int(rnd.Rannyu(1,n-1));
	//cout << r1 <<"   " << r2 << endl;
		for(int i=0; i<n; i++){
		x_new[i]=x[i];
		y_new[i]=y[i];
	
	}
	appo=x_new[r1];
	x_new[r1]=x[r2];
	x_new[r2]=appo;

	
	appo=y_new[r1];
	y_new[r1]=y_new[r2];
	y_new[r2]=appo;
	

}

void MultiplePermutation(){
	int m = int(rnd.Rannyu(1,(n)));
	int r1,r2;
	double appo;
		for(int i=0; i<n; i++){
		x_new[i]=x[i];
		y_new[i]=y[i];
	
	}
	for(int i=0;i<m;i++){

		r1 = int(rnd.Rannyu(1,n-1));
		r2 = int(rnd.Rannyu(1,n-1));
		//cout << r1 <<"   " << r2 << endl;

		appo=x_new[r1];
		x_new[r1]=x_new[r2];
		x_new[r2]=appo;

		
		appo=y_new[r1];
		y_new[r1]=y_new[r2];
		y_new[r2]=appo;
	}
}

void ContinuosPermutation(){
	int s,m,s2;
	s = int(rnd.Rannyu(1,n));
	s2 = int(rnd.Rannyu(1,n));
	m = int(rnd.Rannyu(1,(n/2)));
	m=1;
	double appo,appo2,appo3,appo4;
	for(int i=0; i<n; i++){
		x_new[i]=x[i];
		y_new[i]=y[i];
	
	}

	for(int i=0; i<m; i++){
		if(((s+i)%n==0) or ((s2+i)%n==0)){
			}
			else{
			appo=x_new[(s+i)%n];
			appo2=x_new[(s2+i)%n];
			appo3=y_new[(s+i)%n];
			appo4=y_new[(s2+i)%n];
			x_new[(s2+i)%n] = appo;
			x_new[(s+i)%n]=appo2;
			y_new[(s2+i)%n] = appo3;
			y_new[(s+i)%n] = appo4;
			}
		}
	

}

void Inversion(){
	int s,m,s2;
	s = int(rnd.Rannyu(1,n));
	m = int(rnd.Rannyu(1,(n/2)));
	double appo[n],appo2[n];
	
	for(int i=0; i<n; i++){
		x_new[i]=x[i];
		y_new[i]=y[i];
		appo[i]=x_new[i];
		appo2[i]=y_new[i];
	}

	for(int i=0; i<m; i++){
		if(((s+i)%n==0) or ((s+m-1-i)%n==0)){
			}
			else{
			x_new[(s+i)%n] = appo[(s+m-1-i)%n];
			y_new[(s+i)%n] = appo2[(s+m-1-i)%n];
			}
		}


} 

void Print(string c){	//stampiamo le posizioni delle città

	ofstream out;
	out.open(c);
	for(int i=0; i < n; i++){
		out << x[i] << setw(12) << y[i] << endl;

		}
	
	out.close();

}















