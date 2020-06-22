#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "random.h"
#include "esercitazione91.h"
#include <fstream>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
using namespace std;

int main(){
	Input(); //legge file da settings.C
	clean(); // pulisce i file di output
	BuildingConfig(); //cerchio o quadrato

	for(int d=0; d<depth; d++){ //evoluzione per d=400 generazioni
		if(d%1==0){cout << "Generation: " << d << endl<<endl;}
		for(int p0=0; p0<pop; p0++){
			
			Mutation(p0);	
		}

		bestl2(); //trova il miglior cammino della generazione
		Crossover(); //fa un crossover e inizia una nuova generazione

	}
	Print("best"+config+".dat"); //stampa i migliori cammini
	

return 0;
}
void bestl2(){
	double a[pop];

	for(int p0=0; p0<pop; p0++){ //trova il min
		a[p0]=Evaluation(p0);
		if (a[p0]<minl2) minl2=a[p0];
		}
	ofstream out;
	out.open("l2"+config+".dat",ios::app); //stampa
	out << minl2 << endl;
	out.close();
	sort(a,a+pop); //ordino il vettore

	double sum=0;
	out.open("l2_besthalf"+config+".dat",ios::app);
	for(int p0=0; p0<pop/2; p0++){
		sum += a[p0];
		}

	out << 2*sum/pop << endl; //media dei primi metà cammini miglori
	out.close();

	//per permutazione equa (vedere multiple permutation)
	best_m=a[0];

}


void clean(){
	ofstream out;
	out.open("l2"+config+".dat");
		out <<"";
	out.close();
	out.open("l2_besthalf"+config+".dat");
		out <<"";
	out.close();


}
void Input(){
	ifstream read;

	//Read seed for random numbers
	rnd.LoadSeed();



	read.open("settings.C");
	cout << endl;
	cout << "***********************************************************************" << endl;
	cout << "You are using a genetic algorithm to solve a Traveling Salesman Problem " << endl<<endl;

	read >> config;
	cout << "The cities configuration is a " << config << endl;

	//read >> n;
	cout << "There are "<< n << " cities" << endl;

	cout << "There are " << pop <<" vector of " << n << " cities" << endl;

	read >> depth;
	cout << "The depth of the evolving algorithm is " << depth << " steps"<<endl;
	read >> pm;
	read >> pc;
	cout << "Probability of crossover= " << pc << endl;
	cout << "Probability of mutation= " << pm << endl;
}

void BuildingConfig(){

	if(config == "circle") {
		for(int i=0; i<n; i++){
			x[0][i] = rnd.Rannyu(-1,1);
			y[0][i] = pow((-1),i)*sqrt(1-x[0][i]*x[0][i]);
			
 			}
		}
	if(config == "square") {
		for(int i=0; i<n; i++){
			x[0][i] = rnd.Rannyu(0,1);
			y[0][i] = rnd.Rannyu(0,1);
			}  
		}


	vector <int> vector; //uso vector per mescolare i vettori con random_shuffle
	for(int i=0; i<n; i++){
			vector.push_back(i);
			}  

	for (int j=1; j<pop; j++){ //creo una popolazione di vettori mescolati rispetto a quelli iniziali, senza mischiare le x con le y, quindi mescolo degli indici non direttamente x e y
		random_shuffle(vector.begin(),vector.end());
		
		for(int i=0; i<n; i++){

			
			x[j][i]=x[0][vector[i]];
			y[j][i]=y[0][vector[i]];
			
			} 
		
		
		} 



}

double Evaluation(int actual){ //l2 in realtà non l1 come il nome della variabile 
	l1=0;
	for(int i=0; i<n-1;i++){
		l1 += abs(sqrt(pow((x[actual][i]-x[actual][i+1]),2)+pow((y[actual][i]-y[actual][i+1]),2)));
		}
	l1 += abs(sqrt(pow((x[actual][n-1]-x[actual][0]),2)+pow((y[actual][n-1]-y[actual][0]),2)));	
	return l1;
}






void Crossover(){ 
	double probCross,b,b1,c1,c,appo;
	int N,N1,N2;
	
	
	
	
	for(int e=0; e<pop; e=e+2){
		N=int(rnd.Rannyu(1,n));
		
		b1=9999;
		c1=9999;
		probCross=rnd.Rannyu(0.7,1);
		
		//SELECTION
		for(int i=0; i<pop; i++){
			b=rnd.Rannyu(0.9,1)*Evaluation(i);	//altre possibilità? Si, ma così sembra sensato. (Notebook...)
			

			if(b1>b){N1=i;
				b1=b;}
			}
		

		for(int i=0; i<pop; i++){
			
			if(i != N1){
				c=rnd.Rannyu(0,1)*Evaluation(i);
				if(c1>c){
					N2=i;
					c1=c;}
				}
			}
		

	
		
		//CROSSOVER
		if(probCross<pc){
			for(int no=0; no<N; no++){	
				
				x_new[e][no]=x[N1][no];
				x_new[e+1][no]= x[N2][no];

				y_new[e][no]=y[N1][no];
				y_new[e+1][no]=y[N2][no];
				}
			int appo;
			int k=0;
			for(int i=0; i<n; i++){
				appo=0;
				for(int no=0; no<N; no++){
					if(x[N2][i]==x_new[e][no]){appo=1;}
					}
				if(appo==0){x_new[e][N+k]=x[N2][i];
					    y_new[e][N+k]=y[N2][i];
						k++;}
				
					
				}
			k=0;
			for(int i=0; i<n; i++){
				appo=0;
				for(int no=0; no<N; no++){
					if(x[N1][i]==x_new[e+1][no]){appo=1;}
					}
				if(appo==0){x_new[e+1][N+k]=x[N1][i];
					    y_new[e+1][N+k]=y[N1][i];
						k++;}
				
					
				}
				
			
			}
		else {
			
			for(int no=0; no<n; no++){
				
				x_new[e][no]=x[N1][no];
				x_new[e+1][no]=x[N2][no];
				y_new[e][no]=y[N1][no];
				y_new[e+1][no]=y[N2][no];
				}
		
		}


	}
	Generation(); //funzione per aggiornare la generazione
	
}









void Generation(){//funzione per aggiornare la generazione

	for(int j=0; j<pop; j++){
		
		for(int no=0; no<n; no++){
				
				x[j][no]=x_new[j][no];
				y[j][no]=y_new[j][no];
				}
		}
}




void Mutation(int actual){ // da mantenere fermo il primo elemento , funzione generale di mutazione

	double a,r;
	a=rnd.Rannyu(0,1);
	if (a < pm){ r=int(rnd.Rannyu()*4);}
	
	if(r == 0){PairPermutation(actual);}
	if(r == 1){MultiplePermutation(actual);}
	if(r == 2){ContinuosPermutation(actual);}
	if(r == 3){Inversion(actual);}
}

void PairPermutation(int actual){ //primo elemento non permuta, sceglie due punti di un vettore di città  e li permuta

	int r1,r2;
	double appo;
	r1 = int(rnd.Rannyu(1,n-1));
	r2 = int(rnd.Rannyu(1,n-1));
	//cout << r1 <<"   " << r2 << endl;

	appo=x[actual][r1];
	x[actual][r1]=x[actual][r2];
	x[actual][r2]=appo;

	
	appo=y[actual][r1];
	y[actual][r1]=y[actual][r2];
	y[actual][r2]=appo;
	

}

void MultiplePermutation(int actual){//simile a quella di prima, ma ne sceglie m sulla base di quanto un cammino ne ha bisogno :) (alto se la sua L è alta, basso se la sua L è bassa, rispetto al miglior cammino)
	double f= Evaluation(actual)/best_m;
	
	int m = int(rnd.Rannyu(1,(n*f/8)));
	//cout << m << endl;
	int r1,r2; //città che permutano
	double appo;
	for(int i=0;i<m;i++){

		r1 = int(rnd.Rannyu(1,n-1)); //toglie il primo elemento dalla mutazione
		r2 = int(rnd.Rannyu(1,n-1)); //toglie il primo elemento dalla mutazione
		

		appo=x[actual][r1];
		x[actual][r1]=x[actual][r2];
		x[actual][r2]=appo;

		
		appo=y[actual][r1];
		y[actual][r1]=y[actual][r2];
		y[actual][r2]=appo;
	}
}

void ContinuosPermutation(int actual){ //permuta m città vicine con altre m città vicine a partire da s e s2
	int s,m,s2;
	s = int(rnd.Rannyu(1,n));
	s2 = int(rnd.Rannyu(1,n));
	m = int(rnd.Rannyu(1,(n/2)));
	m=1;
	double appo,appo2,appo3,appo4;


	for(int i=0; i<m; i++){
		if(((s+i)%n==0) or ((s2+i)%n==0)){//toglie il primo elemento dalla mutazione
			}
			else{
			appo=x[actual][(s+i)%n];
			appo2=x[actual][(s2+i)%n];
			appo3=y[actual][(s+i)%n];
			appo4=y[actual][(s2+i)%n];
			x[actual][(s2+i)%n] = appo;
			x[actual][(s+i)%n]=appo2;
			y[actual][(s2+i)%n] = appo3;
			y[actual][(s+i)%n] = appo4;
			}
		}
	

}

void Inversion(int actual){//inverte m città a partire da s
	int s,m;
	s = int(rnd.Rannyu(1,n));
	m = int(rnd.Rannyu(1,(n/2)));
	double appo[n],appo2[n];
	for(int i=0; i<n; i++){
		appo[i]=x[actual][i];
		appo2[i]=y[actual][i];
	}

	for(int i=0; i<m; i++){
		if(((s+i)%n==0) or ((s+m-1-i)%n==0)){//toglie il primo elemento dalla mutazione
			}
			else{
			x[actual][(s+i)%n] = appo[(s+m-1-i)%n]; //modulo n necessario
			y[actual][(s+i)%n] = appo2[(s+m-1-i)%n];
			}
		}


} 

void Print(string c){	
	double b=9999;  
	int best;
	for(int i=0; i<pop; i++){
		if(b>Evaluation(i)){
			b=Evaluation(i);
			best=i;
			
			}
		if(i%10==0){cout << " Actual best evaluation: " << b << endl;}
		}
	cout << "The super best bestissimo è........: " << best<<endl;
	
	ofstream out;
	out.open(c);

	for(int i=0; i < n; i++){
		out << x[best][i] << setw(12) << y[best][i] << endl;
		}
	
	out.close();

}















