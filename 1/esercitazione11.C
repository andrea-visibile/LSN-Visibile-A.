#include<iostream>
#include <fstream>
#include<cmath>
#include "random.h"
#include "esercitazione11.h"
using namespace std;
/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Esercizio 01.1
/////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	
	//seed 
	rnd.LoadSeed();
	
	//per cambiare il numero di tiri o di blocchi usare esercitazione11.h

	///////////////////////////////// Punto a

	//estraggo i numeri casuali
	for (int i=0; i<M; i++){
		throws[i]=rnd.Rannyu();
		}


	//calcolo le medie negli N blocchi (ovvero ottengo N "dati")
	for (int j=0; j<N; j++){ 
		double sum = 0; //variable somma a zero
		for (int k=0; k<L;k++){ // calcolo gli N dati con L numeri per ognuno
			sum += throws[j*L+k];
			}
		av1[j] = sum/L; //faccio la media qui
		av2[j] = pow(av1[j],2); //sto facendo il quadrato del mio "dato"
		}
	for (int i=0; i<N; i++){ //a zero i vettori in cui accumulo
		cum1[i]=0;
		cum2[i]=0;
	}

	//calcolo la media al variare del numero di blocchi
	for (int i=0; i<N; i++){ 
		for (int j=0; j<i+1; j++){ 
			cum1[i] += av1[j]; //sommo i+1 medie (voglio sommarne almeno una, zero non ha senso; i è zero invece perchè è l'indice del v
			cum2[i] += av2[j]; //sommo i+1 quadrati delle medie
			}
		cum1[i] /= (i+1); //faccio la media di i+1 blocchi,vettore in cui c'è l'andamento del valor medio al variare del numero di blocchi 
		cum2[i] /= (i+1); //uguale ma al quadrato
		err[i]= error( cum1, cum2, i); //vettore di errori sulle medie
		}


	//salvo su file medie.dat
	ofstream fileout; 
	fileout.open("medie.dat");
	for (int i=0; i<N; i++){ 
		fileout << cum1[i] <<endl; //salvo le medie
		}
	for (int i=0; i<N; i++){ 
		fileout << err[i] <<endl; //salvo gli errori
		}
	fileout.close();
	

/////////////////////////////**************************************************************/////////////////////////


	/////////////////////// Punto b

	//calcolo le varianze negli N blocchi (ovvero ottengo N "dati")
	for (int j=0; j<N; j++){ 
		double sum = 0;
		for (int k=0; k<L;k++){
			sum += pow(throws[j*L+k]-0.5,2); //uso la media del blocco o 0.5? 0.5!
			}
		av3[j] = sum/L; //media (delle varianze nel blocco)=dato
		av4[j] = pow(av3[j],2); //sto facendo il quadrato del mio "dato"
		}
	for (int i=0; i<N; i++){
		cumvar1[i]=0;
		cumvar2[i]=0;
	}
	 //calcolo la media al variare del numero di blocchi
	for (int i=0; i<N; i++){
		for (int j=0; j<i+1; j++){ 
			cumvar1[i] += av3[j]; //sommo i+1 varianze
			cumvar2[i] += av4[j]; //sommo i+1 quadrati delle varianze
			}
		cumvar1[i] /= (i+1); //faccio la media di i+1 blocchi
		cumvar2[i] /= (i+1); 
		errvar[i]= error( cumvar1, cumvar2, i); //vettore di errori sulle varianze
		}

	//salvo su file varianze.dat
	fileout.open("varianze.dat");
	for (int i=0; i<N; i++){ 
		fileout << cumvar1[i] <<endl; //salvo le varianze
		}
	for (int i=0; i<N; i++){ 
		fileout << errvar[i] <<endl; //salvo gli errori
		}
	fileout.close();


/////////////////////////////**************************************************************/////////////////////////

	/////////////////// Punto c

	
	
	for(int p=0; p<J; p++){  //devo calcolarmi cento chi quadri
		//inizializzo a zero il vettore
		for (int i=0; i<J; i++){ //uso sempre J, ma in realtà sarebbe M (sono entrambi a 100)
			counter[i]=0;
			}
		//calcolo quanti i elementi degli n estratti finiscono nel k-esimo intervallo
		for (int i=0; i<n; i++){  //devo estrarre mille numeri
			throws[i]=rnd.Rannyu();
			for(int k=0; k<J;k++){
				if (throws[i]>=(k*1.0/J) and throws[i]<((k+1)*1.0/J)){
					counter[k] +=1;				
					}
				}
			}
		chi[p]=0;	
		for(int j=0; j<J; j++) {  //alla fine calcolo chi
			chi[p] += (pow(counter[j]-(n/J),2))*1.0/(n/J);
			}
		}
	
	
	//salvo i chi quadri su "chi.dat"
	fileout.open("chi.dat");
	for (int i=0; i<J; i++){ 
		fileout << chi[i] << endl; 
		}
	fileout.close();
	
	//file di dati:
	// medie.dat ha 100 medie e 100 errori
	// varianze.dat idem con varianze
	// chi.dat 100 chi quadri
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
		return 0; //la prima deve essere zero
		}
	else {
		err=av2[n]-pow(av1[n],2);
		return sqrt(err/n);
		}
	}


//AV



