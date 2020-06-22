#include<iostream>
#include <fstream>
#include<cmath>
#include "random.h"
#include "esercitazione12.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////
                                     //Esercizio 01.2
/////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	
	rnd.LoadSeed(); //funzione implementata in random
	
	//tre file diversi per i tre tipi di dadi
	ofstream out1,out2,out3;
	out1.open("unif.dat");
	out2.open("exp.dat");
	out3.open("lorentz.dat");

	for(int j:a){		//somma j numeri casuali, per i quattro j necessari:1,2,10,100)
		for(int i=0; i<n; i++){  //calcola l'i-esimo media con j numeri casuali per ogni media
			x=0;
			e=0;
			l=0;
			for(int f=0; f<j;f++){ 	//f è l'f-esimo dato che formerà l'i-esima media, ne devo estrarre j 
				x += rnd.Rannyu(); //unif
				e += rnd.RandomExp(rnd.Rannyu(),1); //exp
				l += rnd.RandomLorentz(rnd.Rannyu(),0,1); //lorentziana
				}
			throws[i] = x/(j*1.); //unif 
			throwsExp[i] = e/(j*1.) ; //exp
			throwsLor[i] = l/(j*1.); //lor
			//stampo su file
			out1 << throws[i]<<endl; 
			out2 << throwsExp[i]<<endl; 
			out3 << throwsLor[i]<<endl;
			}
 	 
		}
	out3.close();
	out2.close();
	out1.close();

	//ogni file ha n=10000 medie *4 ( i 1,2,10,100 dati per media, in questo ordine)
	rnd.SaveSeed();


	return 0;
	}

/////////////////////////////////////////////////////////////////////////////////////////////
                                     // 
/////////////////////////////////////////////////////////////////////////////////////////////
//AV





