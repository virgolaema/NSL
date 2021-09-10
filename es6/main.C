/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "../lib.h"
#include "../random.h"

using namespace std;

class Ising{
    public:
    Random rnd;
    double * s; //spin of the particles
    int nstep, nblk, metro, nT = 21; //nT is n of simulations at different T
    double stepT;
    double accepted,attempted;
    double * stima_u, * stima_c,* stima_m,* stima_x; 
    int nspin;
    double beta,temp,J,h;
    Ising();
    void Input(void);
    void Reset(int);
    void Averages();
    void Move();
    void ConfFinal(void);
    void Measure(int);
    void Expectations(int);
    double Boltzmann(int, int);
    int Pbc(int);
    void setTemp (double t){temp = t;}
};

int main(){ 
  Ising sim;
  sim.Input(); //Inizialization
  for (int k = 0; k < sim.nT; k++){
    cout << "Computing for T = " << sim.temp << endl; 
    for(int iblk = 0; iblk < sim.nblk; ++iblk){
      sim.Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= sim.nstep; ++istep){
        sim.Move();
        sim.Measure(iblk);
      }
      sim.Expectations (iblk);
    }
    sim.Averages();   //Print results for current block
  }

  sim.ConfFinal(); //Write final configuration
  return 0;
}

Ising :: Ising (){
  Random appo;
  rnd = appo;
  rnd.Initialize(rnd);
}

void Ising::Input(void){
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  cout << "Initial temperature = " << temp << endl;
  beta = 1./temp;
  double tempf;
  ReadInput >> tempf;
  cout << "Final temperature = " << temp << endl;
  cout << "Steps between in and fin: " << nT << endl;

  stepT =   (tempf - temp) / ((double) nT-1);

  cout << "Step of variation of T: " << stepT << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;
  s = new double [nspin]; 

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

  int M = 0;
  ReadInput >> M;
  nstep = M/nblk;

  if(metro==1) cout << "The program performs Metropolis moves" << endl;
  else cout << "The program performs Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  int prec = 0;
  ReadInput >> prec;
  ReadInput.close();

//initial configuration, choose randomly sign of the spin
  if (prec == 0){
    for (int i=0; i < nspin; ++i)  {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
    

  if (prec == 1){
    cout << "Using old configuration of spins" << endl; 
    ifstream confprec;
    confprec.open("config.0");
    for (int i=0; i < nspin; i++) {
      confprec >> s[i];
    }
    confprec.close();
  }
  
  cout << "Initial config: ";
  for(int f = 0; f < nspin; f++) cout << s[f] << " ";
  cout << endl;
}

void Ising::Move(){
  beta = 1./temp;
  int spin_try;
  for(int i=0; i < nspin; ++i) { // try flipping nspin times, but may be more/less
  //Select randomly a particle 
  spin_try = (int)(rnd.Rannyu()*nspin);
  attempted++; 

  if(metro==1){ //Metropolis
    double uold = 0; //old energy of the configuration
    for (int k=0; k<nspin; ++k) uold += -J * s[k] * s[Pbc(k+1)] - h * s[k];
    
    double unew = 0; //new energy of the configuration
    s[spin_try] *= -1; //spin flip, just an attempt
    for (int k=0; k<nspin; ++k) unew += -J * s[k] * s[Pbc(k+1)] - h * s[k];

    double A = min(1., exp (-beta * (unew - uold)));
    if (rnd.Rannyu() > A) s[spin_try] *= -1; //reject the flip, back to original value
    else {accepted++;} //already flipped, just ++ the count
  }

  else if (metro == 0){ //Gibbs sampling
    double sum_spin = s[Pbc(spin_try-1)] + s[Pbc(spin_try+1)];
    double E_gibbs = - (2. * J * sum_spin + 2.*h);
    double gibbs_p = 1./(1. + exp(beta * E_gibbs));
    //cout << "probabilità gibbs  " << gibbs_p << endl;

    if (gibbs_p > rnd.Rannyu() ) s[spin_try] = 1; 
    else {s[spin_try] = -1;} 
    accepted++; //gibbs accepts 100% 
    }
  }
}

void Ising::Measure(int iblk){
  double u[nspin], appo_x=0, appo_u=0;
  for (int i=0; i<nspin; ++i)  { 
    u[i] = -J * s[i] * s[Pbc(i+1)] -  h * s[i]; 
    appo_u += u[i];
    appo_x += s[i];
    stima_m[iblk] += s[i];
  }
  stima_u [iblk] += appo_u; 
  stima_c[iblk] += pow (appo_u,2);
  stima_x[iblk] += pow (appo_x,2);
}

void Ising::Expectations(int iblk){ //computing exp values for quantities
  stima_u [iblk] /= nstep;
  stima_c [iblk] = pow (beta,2) *(stima_c[iblk] / nstep - pow(stima_u [iblk],2));
  stima_x [iblk] = stima_x [iblk] * beta / nstep;
  stima_m [iblk] /= nstep;  
}

void Ising::Averages(){
    
  ofstream Ene, Heat, Mag, Chi;
    
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
  cout << "----------------------------" << endl << endl;
  string appo;
  if (metro == 1) appo = "M";
  else if (metro == 0){ appo = "G";}

  if (h==0){
    Ene.open("u"+appo+".out",ios::app);  
    Ene << temp << "," << mean(nblk, stima_u) << "," << sqrt(variance_blocks(nblk, stima_u))  << endl;
    Ene.close();
    
    Heat.open("c"+appo+".out",ios::app);  
    Heat << temp << "," << mean(nblk, stima_c) << "," << sqrt(variance_blocks(nblk, stima_c))<< endl;
    Heat.close();

    Chi.open("x"+appo+".out",ios::app);  
    Chi << temp << "," << mean(nblk, stima_x) << "," << sqrt(variance_blocks(nblk, stima_x)) << endl;
    Chi.close();
  }

  if(h!=0){
    Mag.open("m"+appo+".out",ios::app);
    Mag << temp << "," << mean(nblk, stima_m) << "," << sqrt(variance_blocks(nblk, stima_m)) << endl;
    Mag.close();
  }

 //change temp and reset to 0 arrays
  temp += stepT;   
  for (int j = 0; j < nblk; j++){
   stima_u [j] = 0;
   stima_c [j] = 0;
   stima_x [j] = 0;
   stima_m [j] = 0;
 }
}

void Ising::ConfFinal(void){
  ofstream WriteConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i) WriteConf << s[i] << endl;
  WriteConf.close();
  cout << "config.final written" << endl;
  rnd.SaveSeed();
}

void Ising::Reset(int iblk){
  if (iblk%5 == 0) cout << "Step " << iblk << endl;
  attempted = 0;
  accepted = 0;
}

int Ising::Pbc(int i){
    if(i >= nspin) i -= nspin;
    else if(i < 0) i += nspin;
    return i;
}

double Ising::Boltzmann(int sm, int ip){return -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
