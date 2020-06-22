#include<iostream>
#include<fstream>
#include<cmath>
#include "random.h"
#include "esercitazione21.h"
using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Esercizio 02.1
/////////////////////////////////////////////////////////////////////////////////////////////

int main(){

	rnd.LoadSeed();

	//estraggo i numeri casuali
	for (int i=0; i<M; i++){
		x = rnd.Rannyu(); //estratto unif
		y = (-2 + sqrt(4-4*x))/(-2); //cumulativa inv
	
		unif[i] = (M_PI/2)*cos(M_PI*x/2); //uniform sampling
		imp[i] = (M_PI/2)*cos(M_PI*y/2)/(2*(1-y)); //importance sampling
		
		}


	//N blocchi
	for (int j=0; j<N; j++){ 
		double sum = 0;
		double sum2 = 0;
		for (int k=0; k<L;k++){ // calcolo gli N dati con L numeri per ognuno
			sum += unif[j*L+k]; //unif
			sum2 += imp[j*L+k]; //imp
			}
		av1[j] = sum/L; //unif
		av2[j] = pow(av1[j],2); //unif quadro
		av3[j] = sum2/L; //imp
		av4[j] = pow(av3[j],2); //imp quadro
		}

	//metto a zero i vettori
	for (int i=0; i<N; i++){ 
		cum1[i] = 0;
		cum2[i] = 0;
		cum3[i] = 0;
		cum4[i] = 0;

		}


	//calcolo la media al variare del numero di blocchi
	for (int i=0; i<N; i++){ 
		for (int j=0; j<i+1; j++){ 
			cum1[i] += av1[j]; 
			cum2[i] += av2[j]; 
			cum3[i] += av3[j]; 
			cum4[i] += av4[j]; 	
			
			}

		cum1[i] /= (i+1);  
		cum2[i] /= (i+1); 
		err[i]= error( cum1, cum2, i); //errore unif
		cum3[i] /= (i+1);  
		cum4[i] /= (i+1); 
		err2[i]= error( cum3, cum4, i);  //errore imp

		}

	//salvo su file medie.dat
	ofstream fileout, fileout2; 
	fileout.open("uniform_integral.dat");
	fileout2.open("integral.dat");
	for (int i=0; i<N; i++){ 
		fileout << cum1[i] <<endl; //unif
		fileout2 << cum3[i] <<endl; //imp
		}
	for (int i=0; i<N; i++){ 
		fileout << err[i] <<endl; //salvo gli errori
		fileout2 << err2[i] <<endl; 
		}
	fileout.close();


	rnd.SaveSeed();
	return 0;
	}


/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Funzioni
/////////////////////////////////////////////////////////////////////////////////////////////

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

//AV
