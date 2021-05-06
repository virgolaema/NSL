/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/*#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>*/
#include "main.h"

using namespace std;

int main(){ 
  rnd.Initialize(rnd);
  Input(); //Inizialization

  for(int iblk=1; iblk <= nblk; ++iblk){
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move(metro);
      Measure(iblk, istep);
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  delete []stima_u,stima_x, stima_m, stima_c;
  return 0;
}


void Input(void){
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1./temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;
  stima_u = new double [nblk]{};
  stima_c = new double [nblk]{};
  stima_x = new double [nblk]{};
  stima_m = new double [nblk]{};

  ReadInput >> nstep;

  if(metro==1) cout << "The program performs Metropolis moves" << endl;
  else cout << "The program performs Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  n_props = 4; //Number of observables

//initial configuration, choose randomly sign of the spin
  for (int i=0; i < nspin; ++i)  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure(1,1); //compute for first block

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << stima_u[1] << endl;
}

void Move(int metro){
  int spin_try;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i < nspin; ++i) { //nspin times, but may be more/less
  //Select randomly a particle 
    spin_try = (int)(rnd.Rannyu()*nspin);
    attempted++; //useful only for metropolis
    if(metro==1){ //Metropolis
      double uold = 0; //old energy of the configuration
      for (int k=0; k<nspin; ++k) uold += -J * s[k] * s[Pbc(k+1)] - 0.5 * h * (s[k] + s[Pbc(k+1)]); //internal energy, obtained by summing i anf i+1, not both sides
      
      double unew = 0; //new energy of the configuration
      s[spin_try] *= -1; //spin flip, just an attempt
      for (int k=0; k<nspin; ++k) unew += -J * s[k] * s[Pbc(k+1)] - 0.5 * h * (s[k] + s[Pbc(k+1)]);

      double A = min(1., metro_p(unew) / metro_p(uold));
      if (rnd.Rannyu() > A) s[spin_try] *= -1; //reject the flip, back to original value
      else {accepted++;} //already flipped, just ++ the count
    }

    else if (metro == 0){ //Gibbs sampling
      double sum_spin = 0;
      for (int k = 0; k < nspin; k++) sum_spin += s[i];
      double gibbs_p = 1./exp(2 * beta * J * sum_spin);

      if (rnd.Rannyu() < gibbs_p) s[spin_try] *= -1; //flip
      accepted++; //gibbs accepts 100% 
    }
  }
}

double Boltzmann(int sm, int ip){return -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;}

void Measure(int iblk, int istep){
  int bin; //what use is this for
  double u[nspin];

  ofstream val_u, val_c, val_x, val_m;
  val_u.open("u_values.0", ios::app);
  val_c.open("c_values.0", ios::app);
  val_x.open("x_values.0", ios::app);
  val_m.open("m_values.0", ios::app);

//cycle over spins
  for (int i=0; i<nspin; ++i)  { //shouldnt 0.5 go in the first term
    u[i] = -J * s[i] * s[Pbc(i+1)] -  h * s[i]; //internal energy, obtained by summing i anf i+1, not both sides
    stima_u[iblk] += u[i];
    stima_m[iblk] += s[i];
    stima_x[iblk] += s[i];
  }
  stima_u[iblk] /= nspin;
  stima_c[iblk] = beta*beta * variance_blocks(nspin, u) * (nspin-1); 
  stima_x[iblk] = stima_x[iblk] * stima_x[iblk] / nspin * beta;
  stima_m[iblk] /= nspin; //h=0.02 how

  if (iblk == 1 & istep <= 200){
    val_u << (iblk-1)*nstep + istep << "," << stima_u[iblk] << endl;
    val_c << (iblk-1)*nstep + istep << "," << stima_c[iblk] << endl;
    val_x << (iblk-1)*nstep + istep << "," << stima_x[iblk] << endl;
    val_m << (iblk-1)*nstep + istep << "," << stima_m[iblk] << endl;
  }
}

void Averages(int iblk){
    
  ofstream Ene, Heat, Mag, Chi;
    
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
  cout << "----------------------------" << endl << endl;
    
  Ene.open("output.ene.0",ios::app);  
  Ene << iblk << "," << mean(iblk, stima_u) << "," << sqrt(variance_blocks(iblk, stima_u)) << "," << endl;
  Ene.close();
  
  Heat.open("output.heat.0",ios::app);  
  Heat << iblk << "," << mean(iblk, stima_c) << "," << sqrt(variance_blocks(iblk, stima_c)) << "," << endl;
  Ene.close();
  
  Mag.open("output.mag.0",ios::app);  
  Mag << iblk << "," << mean(iblk, stima_m) << "," << sqrt(variance_blocks(iblk, stima_m)) << "," << endl;
  Mag.close();
  
  Chi.open("output.chi.0",ios::app);  
  Chi << iblk << "," << mean(iblk, stima_x) << "," << sqrt(variance_blocks(iblk, stima_x)) << "," << endl;
  Chi.close();
}

void ConfFinal(void){
  ofstream WriteConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i) WriteConf << s[i] << endl;
  WriteConf.close();
  rnd.SaveSeed();
}

void Reset(int iblk){
  attempted = 0;
  accepted = 0;
}

double metro_p (double u){
  double prob = exp (-u*beta); //missing normalization
  return prob;
}

int Pbc(int i){
    if(i >= nspin) i -= nspin;
    else if(i < 0) i += nspin;
    return i;
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
