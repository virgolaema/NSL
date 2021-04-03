/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   //di seguito, adatto codice per generare numeri secondo i miei desideri

  
   int nGen = 0;
   int which = 0;
   
   cout << "Con questo programma, puoi generare dei numeri random con varie distribuzioni." << endl:
   cout << "1 - uniformemente tra 0 e 1" << endl;
   cout << "2 - secondo una distribuzione gaussiana" << endl;

   cout << "Indica ora quale distribuzione ti interessa, inviando il numero corrispondente:" << endl;
   cin >> which;
   
   cout << "Generare quanti dati?" << endl; 
   cin >> nGen;

   //generatore unif

   if (which == 1){
      
      ofstream output("unif.txt");
      for(int i=0; i<nGen; i++){
         output << rnd.Rannyu() << endl;
      }
      output.close();
   }

   //generatore gaus
   if (which == 2){

   }


   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
