/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie,iw;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
void BlockStat(std::string filein, std::string fileout);
void block_error(double* error, double* sum_prog, double* x1, double* x2, int n);
double Error(double, double, int);

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
//g(r)
const int nbins=100;
double bin_size;
double hist_g[nbins];
double hist_ave[nbins];
double stima_g;
double err_gdir[100];
double glob_av[100];
double glob_av2[100];
int nblk=100;
double norm=0;
// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

//average
double blk_av[2];
double global_m[2];
double global_m2[2];
double walker_av[2];


// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Reset(int );
void Accumulate();
void Average(int );
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
