#include<iostream>
#include<fstream>
#include<cmath>
#include "random.h"
#include "esercitazione31.h"
#include <iomanip>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Esercizio 03.1
/////////////////////////////////////////////////////////////////////////////////////////////



int main(){
	rnd.LoadSeed();

	for(int i=0; i<M; i++){
		S = S_0 * exp((r - 0.5*pow(vol,2))*T+vol*rnd.Gauss(0,T));
		throws[i] = max(0., S-K)*exp(-r*T);
		throws2[i] = max(0., K-S)*exp(-r*T);
		}	

	//calcolo le medie negli N blocchi (ovvero ottengo N "dati")
	for (int j=0; j<N; j++){ 
		double sum = 0;
		double sum2 =0;
		for (int k=0; k<L;k++){ // calcolo gli N dati con L numeri per ognuno
			sum += throws[j*L+k];
			sum2 += throws2[j*L+k];	
			}
		av1[j] = sum/L;
		av2[j] = pow(av1[j],2); //sto facendo il quadrato del mio "dato"
		av3[j] = sum2/L;
		av4[j] = pow(av3[j],2); 
		}
	
	//calcolo la media al variare del numero di blocchi
	for (int i=0; i<N; i++){ 
		for (int j=0; j<i+1; j++){ 
			cum1[i] += av1[j]; //sommo i+1 medie
			cum2[i] += av2[j]; //sommo i+1 quadrati delle medie
			cum3[i] += av3[j]; //sommo i+1 medie
			cum4[i] += av4[j]; //sommo i+1 quadrati delle medie
			}
		cum1[i] /= (i+1); //faccio la media di i+1 blocchi,vettore in cui c'è l'andamento del valor medio al variare del numero di blocchi 
		cum2[i] /= (i+1); 
		cum3[i] /= (i+1); //faccio la media di i+1 blocchi,vettore in cui c'è l'andamento del valor medio al variare del numero di blocchi 
		cum4[i] /= (i+1); 		
		err[i]= error( cum1, cum2, i); //vettore di errori sulle medie
		err2[i]= error( cum3, cum4, i); //vettore di errori sulle medie
		}

	//salvo su file medie.dat
	ofstream fileout; 
	fileout.open("direct_call.dat");
	for (int i=0; i<N; i++){ 
		fileout << cum1[i]<< setw(12)<< err[i] <<endl; //salvo le medie
		}

	fileout.close();
	
	fileout.open("direct_put.dat");
	for (int i=0; i<N; i++){ 
		fileout << cum3[i] <<setw(12)<< err2[i] <<endl; //salvo le medie
		}

	fileout.close();



	

	/////////////////////////////////////************************************////////////////////////////////////////
	//punto b 
	double s_t[100];
	
	s_t[0] = S_0;
	for(int i=0; i<M; i++){
		for(int j=1; j<100; j++){ //DeltaT=1/j
			s_t[j] = s_t[j-1]*exp((r-0.5*pow(vol,2))*(1/100.)+vol*rnd.Gauss(0,1)*sqrt(1/100.)); //cambiare 100 se si vuole generalizzare
			
			}
		S = s_t[99]; // il prezzo finale
		throws[i] = max(0., S-K)*exp(-r*T); //prezzo dell'opzione 
		throws2[i] = max(0., K-S)*exp(-r*T);
		
		}	

	//calcolo le medie negli N blocchi (ovvero ottengo N "dati")
	for (int j=0; j<N; j++){ 
		double sum = 0;
		double sum2 = 0;
		for (int k=0; k<L; k++){ // calcolo gli N dati con L numeri per ognuno
			sum += throws[j*L+k];
			sum2 += throws2[j*L+k];	
			}
		av1[j] = sum/L;
		av2[j] = pow(av1[j],2); //sto facendo il quadrato del mio "dato"
		av3[j] = sum2/L;
		av4[j] = pow(av3[j],2); 
		}
	
	//calcolo la media al variare del numero di blocchi
	for (int i=0; i<N; i++){ 
		cum1[i] = 0;
		cum2[i] = 0;
		cum3[i] = 0;
		cum4[i] = 0;
		for (int j=0; j<i+1; j++){ 
			cum1[i] += av1[j]; //sommo i+1 medie
			cum2[i] += av2[j]; //sommo i+1 quadrati delle medie
			cum3[i] += av3[j]; //sommo i+1 medie
			cum4[i] += av4[j]; //sommo i+1 quadrati delle medie
			}
		cum1[i] /= (i+1); //faccio la media di i+1 blocchi,vettore in cui c'è l'andamento del valor medio al variare del numero di blocchi 
		cum2[i] /= (i+1); 
		cum3[i] /= (i+1); //faccio la media di i+1 blocchi,vettore in cui c'è l'andamento del valor medio al variare del numero di blocchi 
		cum4[i] /= (i+1); 		
		err[i]= error( cum1, cum2, i); //vettore di errori sulle medie
		err2[i]= error( cum3, cum4, i); //vettore di errori sulle medie
		}

	//salvo su file medie.dat

	fileout.open("discrete_call.dat");
	for (int i=0; i<N; i++){ 
		fileout << cum1[i]<<setw(12) << err[i] <<endl; //salvo le medie
		}

	fileout.close();
	
	fileout.open("discrete_put.dat");
	for (int i=0; i<N; i++){ 
		fileout << cum3[i] <<setw(12)<< err2[i] <<endl; //salvo le medie
		}
	fileout.close();

	
	rnd.SaveSeed();


	return 0;	
}


/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Funzioni
/////////////////////////////////////////////////////////////////////////////////////////////

// funzione per calcolare l'errore con il metodo a blocchi
double error(double av1[], double av2[], int n) {
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


