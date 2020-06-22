#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include <string>
#include "esercitazione51.h"
#include <iomanip>
using namespace std;

/////////////////////////////////*********************************************//////////////////////////////////////////////
                                              //Esercizio5.1
/////////////////////////////////*********************************************//////////////////////////////////////////////


int main(){
	rnd.LoadSeed();
	settings(); //Fondamentale! FILE DI INPUT!

	for (int i=0; i<N; i++){ //vettori a zero per media a blocchi
		cum1[i]=0;
		cum2[i]=0;
	}

	//calcolo le medie negli N blocchi (ovvero ottengo N "dati")
	for (int j=0; j<N; j++){ 
		double sum = 0;
		for (int k=0; k<L;k++){ // calcolo gli N dati con L numeri per ognuno
			sum += metropolis(x);//muovo con il metropolis
			}
		pos[j][0]=x[0]; //per config 3D
		pos[j][1]=x[1];
		pos[j][2]=x[2];
		av1[j] = sum/L;
		av2[j] = pow(av1[j],2); //sto facendo il quadrato del mio "dato"
		}
	
	//calcolo la media al variare del numero di blocchi
	for (int i=0; i<N; i++){ 
		for (int j=0; j<i+1; j++){ 
			cum1[i] += av1[j]; //sommo i+1 medie
			cum2[i] += av2[j]; //sommo i+1 quadrati delle medie
			}
		cum1[i] /= (i+1); //faccio la media di i+1 blocchi,vettore in cui c'è l'andamento del valor medio al variare del numero di blocchi 
		cum2[i] /= (i+1); 
		err[i]= error( cum1, cum2, i); //vettore di errori sulle medie
		}

	//stampo automaticamente sui file giusti da settings.C
	ofstream out, out2, posz;
	out.open("hidrogen" + orbital + function + ".dat");
	out2.open("hidrogen"+ orbital + function + "_err.dat");
	posz.open("hidrogen_pos"+ orbital + function +".dat"); //per la config 3D

	for(int i=0; i<N; i++){
		out << cum1[i]<< endl;
		out2 << err[i] << endl;
		posz << pos[i][0] << setw(14) << pos[i][1] << setw(14) << pos[i][2] << endl;
		}
	rnd.SaveSeed();

	cout << acc*1./M<<endl;
	cout << endl;
	posz.close();
	out2.close();
	out.close();

}

/////////////////////////////////*********************************************//////////////////////////////////////////////
                                              //funzioni 
/////////////////////////////////*********************************************//////////////////////////////////////////////

//calcolare l'errore nella media a blocchi
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

//distribuzione uniforme
double uniform(double x){

	return (x  + (rnd.Rannyu(-delta,delta)));

}

//distribuzione gaussiana
double gaussian(double x){

	return rnd.Gauss(x,delta);

}

//modulo quadro della funzione di probabilità del livello 100
double hidrogen100(double r, double theta){
	double pH;
	pH = exp(-2.*r)/(M_PI);
	return pH;
}

//modulo quadro della funzione di probabilità del livello 210
double hidrogen210(double r, double theta){
	double pH;
	pH = exp(-r)*r*r/(32*M_PI)*pow(cos(theta),2);
	return pH;
}


//Metropolis!
double metropolis(double x[]){
        double y[3];
	double p,R_1,R_2,px,py,r,theta1,theta2;
	for(int i=0; i<3;i++){
		y[i]= t(x[i]);
		
	}

	R_1=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	R_2=sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
	theta1=acos(x[2]/R_1);
	theta2=acos(y[2]/R_2);

	px=f(R_1,theta1);
	py=f(R_2,theta2);
	p=min(1.,py/px);
	

	r=rnd.Rannyu();
	if (r<=p){
		acc++;
		for(int i=0; i<3; i++) {
			x[i]=y[i];
		        
			}
	}
	
	return R_1;

}

//Funzione per evitare di dover cambiare ogni volta tutti i parametri
//settings.C contiene informazioni sull'orbitale e sulla funzione da usare in metropolis
void settings(){
	ifstream in;
	in.open("settings.C");

	in >> orbital;
	in >> function;
	in >> delta;
	in >> x_0;
	in >> y_0;
	in >> z_0;

	x[0]=x_0;
	x[1]=y_0;
	x[2]=z_0;

	if (function == "uniform") t=&uniform;
	if (function =="gaussian") t=&gaussian;

	if (orbital == "100") f=&hidrogen100;
	if (orbital =="210") f=&hidrogen210;

}

//AV
