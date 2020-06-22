#include<iostream>
#include <fstream>
#include<cmath>
#include "random.h"
#include "esercitazione13.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Esercizio 01.3
/////////////////////////////////////////////////////////////////////////////////////////////

int main(){

	//classe random
	Random rnd;
	rnd.LoadSeed();
	
	//per salvare su file
	ofstream out;
	out.open("buffon.dat");

	//calcolo le medie negli N blocchi (ovvero ottengo N "dati")
	for (int j=0; j<100; j++){ 
		double sum = 0;
		b=0;
		for (int k=0; k<L;k++){ //trovo L valori di pi e calcolo una delle N medie
			
			y = abs(rnd.Rannyu(0,d)); //estraggo il punto di arrivo della sbarra
			a = abs(rnd.Rannyu(0,M_PI/2.)); //estraggo un angolo (per togliere pi potrei estrarre la x e fare pitagora
			b += buffon(d,l,a,y); //calcolo se la sbarra tocca o meno le linee 
			
			}
		sum = pi(b,L, d,l); //tiro L dati per un pi
		av1[j] = sum;	
		av2[j] = pow(sum,2); //sto facendo il quadrato del mio "dato"
		}

	for (int i=0; i<N; i++){ 
		cum1[i]=0;
		cum2[i]=0;
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

	//salvo i dati su file
	for (int i=0; i<N; i++){ 
		out << cum1[i] <<endl; //salvo le medie
		}
	for (int i=0; i<N; i++){ 
		out << err[i] <<endl; //salvo gli errori
		}

	out.close();
	


	//dati su file:
	//buffon.dat ha dentro 100 valori di pi e 100 errori in quest'ordine
	rnd.SaveSeed();
	return 0;
	}


/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Funzioni
/////////////////////////////////////////////////////////////////////////////////////////////
//Funzioni
//funzione per calcolare l'errore con il metodo a blocchi
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

//funzione per calcolare pi greco dati un numero di tiri N_throws nell'esperimento di Buffon, il numero di volte N_hit che la sbarra di lunghezza L cade sulle linee distanti d 
double pi(int N_hit, int N_throws,double D, double L){
	double pi=(2.*L*N_throws)/(N_hit*D);

	return pi; 
	}

//Funzione per implementare l'esperimento di Buffon
int buffon(double d, double L, double a, double y){
	double mod,mod_f;

	mod=y-(L/2.)*cos(a);
	//cout << mod << endl;
	mod_f=y+(L/2.)*cos(a);
	//cout << mod_f << endl<<endl;
	if ((mod<0) or (mod_f>d)) return 1;
	else return 0;
	}

//AV





