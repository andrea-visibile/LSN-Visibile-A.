#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "random.h"
#include "esercitazione91.h"
#include <fstream>

using namespace std;

int main(){
	Input();
	BuildingConfig();
	for(int d=0; d<depth; d++){
		for(int p0; p0<pop; p0++){
			Crossover(p0);
			Mutation(p0);	
			e[p0]=Evaluation(p0);
		}
	}
	Print();
return 0;
}

void Input(){
	ifstream read;
	read.open("settings.C");
	cout << "You are using a genetic algorithm to solve a Traveling Salesman Problem " << endl;

	read >> config;
	cout << " The cities configuration is a " << config << endl;

	//read >> n;
	cout << " There are "<< n << " cities" << endl;
	double z,zz;
	read >> z;
	read >> zz;
	//read >> pop;
	cout << "There are " << pop <<" vector of " << n << " cities" << endl;

	read >> depth;
	cout << "The depth of the evolving algorithm is " << depth << " steps";
	read >> pm;
	read >> pc;
	cout << "Probability of crossover= " << pc << endl;
	cout << "Probability of mutation= " << pc << endl;
}

void Buildingconfig(){
	for(int j=0; j<pop; j++){
		if(config == "circle") {
			for(int i=0; i<n; i++){
				x[j][i] = rnd.Rannyu();
				y[j][i] = sqrt(1-x[j][i]*x[j][i]);
	 			}
			}
		if(config == "square") {
			for(int i=0; i<n; i++){
				x[j][i] = rnd.Rannyu();
				y[j][i] = rnd.Rannyu();
				}
			}
		}


}

double Evaluation(int actual){
	for(int i=0; i<n;i++){
		l1 += abs(pow((x[actual][i]-x[actual][i+1]),2)+pow((y[actual][i]-y[actual][i+1]),2));
		}

	return l1;
}






void Crossover(int actual){
	double a,appo;
	int p,N;
	if (a<pc){ 
	p=int(rnd.Rannyu()*pop);
	N=int(rnd.Rannyu()*n);	
		for(int no=n; no>N; no--){

			appo= x[actual][no];
			x[actual][no]=x[p][no];
			x[p][no]= appo;

			appo= y[actual][no];
			y[actual][no]=y[p][no];
			y[p][no]= appo;
			}
	}


}

void Mutation(int actual){

	double a,r;
	if (a < pm){ r=int(rnd.Rannyu()*4);}
	
	if(r == 0){PairPermutation(actual);}
	if(r == 0){Shift(actual);}
	if(r == 0){ContinuosPermutation(actual);}
	if(r == 0){Inversion(actual);}
}

void PairPermutation(int actual){
	int r1,r2;
	double appo;
	int mut;
	mut = int(rnd.Rannyu()*pop);
	r1 = int(rnd.Rannyu()*n);
	r2 = int(rnd.Rannyu()*n);

	appo=x[actual][r1];
	x[actual][r1]=x[mut][r2];
	x[mut][r2]=appo;

	appo=y[actual][r1];
	y[actual][r1]=y[mut][r2];
	y[mut][r2]=appo;

}

void Shift(int actual){
	int s;
	s=int(rnd.Rannyu()*n);
	double appo[n],appo2[n];
	for(int i=0; i<n; i++){
		appo[i]=x[actual][i];
		appo2[i]=y[actual][i];
	}
	 
	for(int i=0; i<n; i++){
		x[actual][i] = appo[(i+s)%n];
		y[actual][i] = appo[(i+s)%n];
		}
}

void ContinuosPermutation(int actual){
	int s,m;
	s = int(rnd.Rannyu()*n);
	m = int(rnd.Rannyu()*n/2.);
	double appo[n],appo2[n];
	for(int i=0; i<n; i++){
		appo[i]=x[actual][i];
		appo2[i]=y[actual][i];
	}

	for(int i=0; i<m; i++){
		x[actual][s+i] = appo[s+1+(i%m)];
		y[actual][s+i] = appo[s+1+(i%m)];
		}

}

void Inversion(int actual){
	int s,m;
	s = int(rnd.Rannyu()*n);
	m = int(rnd.Rannyu()*n/4.)*2;
	double appo[n],appo2[n];
	for(int i=0; i<n; i++){
		appo[i]=x[actual][i];
		appo2[i]=y[actual][i];
	}

	for(int i=0; i<m/2; i++){
		x[actual][s+i] = appo[s+m-i];
		y[actual][s+i] = appo[s+m-i];
		}


}

void Print(){
	fstream out;
	out.open("best.dat");
	for(int i=0; i < n; i++){
		out << best[i] << endl;
		}

	out.close();

}















