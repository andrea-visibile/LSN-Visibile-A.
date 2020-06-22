#include<iostream>
#include<fstream>
#include<cmath>
#include "random.h"
#include "esercitazione22.h"
#include<iomanip>
using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Esercizio 02.2
/////////////////////////////////////////////////////////////////////////////////////////////


int main () {


	rnd.LoadSeed();
	
	//punto a: reticolo cubico

	for(int i=0; i<step; i++){// reitero l'intero procedimento in base a quanti step voglio
		for(int k=0; k<N;k++){//sommo sui blocchi
			for (int j=0; j<n/N; j++){ //faccio un blocco
				x[0] = 0;
				x[1] = 0;
				x[2] = 0;
					for(int f=0; f<i; f++){
						ax =int(rnd.Rannyu()*6); //scelta dell'asse, estrae 0,1,2,3,4,5 es:2
						if(ax%2 == 0){ //2%2=0 quindi avanza di 1 y
							x[ax/2] += 1;
							}
						else{
							x[ax/2] -=1;
							}
						}
				r[j] = pow(x[0],2) + pow(x[1],2) + pow(x[2],2); //calcolo r2
				}

			double sum = 0;
			for (int j=0; j<n/N; j++){
				sum += r[j]; //sommo tutti gli r^2 del blocco
				}
			av1[i][k] = (sum/(n/N)); //trovo r^2 medio per blocco (k blocco, i step)
			
			av2[i][k] = pow(av1[i][k],2);
			}
		
		}
	//punto b: reticolo continuo

	//analogo a prima
	for(int i=0; i<step; i++){
		for(int k=0; k<N;k++){
			for (int j=0; j<n/N; j++){
					
				x2[0] = 0;
				x2[1] = 0;
				x2[2] = 0;
					for(int f=0; f<i; f++){
						theta = acos(1-2*rnd.Rannyu(0,1));
						phi = rnd.Rannyu(0,1)*2*M_PI;
						x2[0] += sin(theta)*cos(phi);
						x2[1] += sin(theta)*sin(phi);
						x2[2] += cos(theta);
						
						}
					
				r[j] = pow(x2[0],2) + pow(x2[1],2) + pow(x2[2],2); //calcolo r2
				
				}
			double sum2 = 0;
			for (int j=0; j<n/N; j++){
				sum2 += r[j];
				}
			av3[i][k] = (sum2/(n/N)); //media di r2 su istep nel blocco k-esimo
			av4[i][k] = pow(av3[i][k],2);
			}
			
		}

	//Salvo tutto su file
	fstream out, out2;
	out.open("markov1");
	out2.open("markov2");

	for(int i=0; i<step; i++){ 
		for (int i=0; i<N; i++){ 
			cum1[i] = 0;
			cum2[i] = 0;
			cum3[i] = 0;
			cum4[i] = 0;
		}

		for(int k=0;k<N;k++){
			cum1[i]+=av1[i][k];//cubico
			cum2[i]+=av2[i][k];//cubico al quadrato
			cum3[i]+=av3[i][k];//continuo
			cum4[i]+=av4[i][k];//continuo al quadrato
			}
		cum1[i] /=double(N); //media delle medie
		cum2[i] /=double(N);
		cum3[i] /=double(N);
		cum4[i] /=double(N);
		err[i]=error(cum1,cum2,i);//errore
		err2[i]=error(cum3,cum4,i);
		
		out  << sqrt(cum1[i]) << setw(12) << err[i] << endl; //salvo su file in maniera più intelligente
		out2 << sqrt(cum3[i]) << setw(12) << err2[i] << endl;
		}
		
	out.close();
	out2.close();
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

