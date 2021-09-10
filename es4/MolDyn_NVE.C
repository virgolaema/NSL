/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "../random.h"
#include "../lib.h"

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
   
     if (istep == nstep-1) ConfFinal("config.00"); //Write config at step t-dt
  }
  ConfFinal("config.final");         //Write final configuration to restart
  Blocking ();

  return 0;
}

///////////////////////////////////////////////////////////////
//Functions

void Input(){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf0, ReadConf00;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  Random rnd;
  rnd.Initialize(rnd);
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf0.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf0 >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box; //multiply for the dimension of the box, in was in LJ units
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf0.close();

  cout << "Do you want to restart reading also the configuration at the time t-dt?" << endl;
  cout << "Type 1 to confirm, 0 to decline" << endl;

  int opt;
  cin >> opt;

  if (opt == 0){
  //Prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rnd.Rannyu() - 0.5; 
      vy[i] = rnd.Rannyu() - 0.5;
      vz[i] = rnd.Rannyu() - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    
    for (int idim = 0; idim < 3; ++idim) sumv [idim] /= (double) npart; //sum of velocities in every direction, per part
    
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){ //renormalize velocities to have null mean, each direction
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]; //total module^2 of v every particle
    }
    double meanv2 = sumv2 / (double)npart; //for each particle

    fs = sqrt(3 * temp / meanv2);   // fs = velocity scale factor, like a normalization so that v^2=3T holds (m=1)
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta); //will be used in Move
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  //if opt is 1, overwrite the old starting coordinate with values from file 00
  if (opt == 1){
    ReadConf00.open("config.00");
    for (int i=0; i<npart; ++i){
      ReadConf0 >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box; //multiply for the dimension of the box, in was in LJ units
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf00.close();

    cout << "The 00 config has been considered for the step t-dt." << endl;
    cout << "Do you want to rescale the coordinates in order to match the desired T?" << endl;
    cout << "Type 1 to confirm, 0 to decline:" << endl;
    int appo;
    cin >> appo;
    if (appo == 1){
      double v2tot = 0;
      for (int i=0; i<npart; ++i){
        vx[i] = (x[i] - xold[i])/delta;
        vy[i] = (y[i] - yold[i])/delta;
        vz[i] = (z[i] - zold[i])/delta;

        v2tot += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]; //total module^2 of v every particle
      }
      double v2tot_new = 3*temp*npart;
      double cf = v2tot_new / v2tot; //find the correction factor

      for (int i=0; i<npart; ++i){ //now rescale
        xold[i] = Pbc(x[i] - sqrt(cf)* vx[i] *delta);
        yold[i] = Pbc(y[i] - sqrt(cf)* vy[i] *delta);
        zold[i] = Pbc(z[i] - sqrt(cf)* vz[i] *delta);
      } //now the temperature is the wanted one

      cout << "Temperature at the step t-dt corrected." << endl<<endl;

    }
  }
  rnd.SaveSeed();
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) ); //f=a in LJ units
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f = 0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){ //radius where F is considered
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); //kin en. total
  
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  
  return;
}

void Blocking (){
  const double kB = 1.38e-23; //j K^-1
  double epsilon = 1, TSI = 1;

  cout << "Do you want to switch to SI units for the ave_output?" << endl;
  cout << "Type 1 to confirm, 0 to decline" << endl;
  int opt;
  cin >> opt;

  if(opt == 1){
    epsilon = 120*kB;
    TSI = 120; 
  }

  int L = 10;
  double val = 0;
  ofstream ave_epot, ave_ekin, ave_etot, ave_temp;
  ave_epot.open ("ave_epot.out");
  ave_ekin.open ("ave_ekin.out");
  ave_etot.open ("ave_etot.out");
  ave_temp.open ("ave_temp.out");

  for (int i = 1; i < nstep/10./L; i++){ 
    ifstream Epot, Ekin, Etot, Temp;
    Epot.open("output_epot.dat");
    Ekin.open("output_ekin.dat");
    Temp.open("output_temp.dat");
    Etot.open("output_etot.dat");
    
    int throws = (i+1)*L; //number of throws at this step
    int Nblocks  = i+1;   //number of blocks at this step
    double  ETOT [Nblocks] = {}, EKIN [Nblocks] = {},EPOT [Nblocks] = {}, TEMP [Nblocks] = {};
    for (int j = 0; j < Nblocks ; j++){
      for (int k = 0; k < L; k++){
        Epot >> val;
        EPOT [j] += val;
        Ekin >> val;
        EKIN [j] += val; 
        Etot >> val;
        ETOT [j] += val;
        Temp >> val;
        TEMP [j] += val;       
      }
      ETOT[j] /= L/epsilon;
      EKIN[j] /= L/epsilon;
      EPOT[j] /= L/epsilon;
      TEMP[j] /= L/TSI;
    }
    ave_epot << Nblocks << "," << mean(Nblocks, EPOT) << "," << variance_blocks (Nblocks, EPOT) << endl;
    ave_ekin << Nblocks << "," << mean(Nblocks, EKIN) << "," << variance_blocks (Nblocks, EKIN) << endl;
    ave_etot << Nblocks << "," << mean(Nblocks, ETOT) << "," << variance_blocks (Nblocks, ETOT) << endl;
    ave_temp << Nblocks << "," << mean(Nblocks, TEMP) << "," << variance_blocks (Nblocks, TEMP) << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
  }
  
  ave_epot.close();
  ave_ekin.close();
  ave_etot.close();
  ave_temp.close();
  return;
}

void ConfFinal(std::string filename){ //Write final configuration
  ofstream WriteConf;
  WriteConf.open(filename);
  for (int i=0; i<npart; ++i)   WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;
  string number = to_string(nconf);
  if (nconf < 10) number = to_string(0) + number;
  if (nconf < 100) number = to_string(0) + number;
  if (nconf < 1000) number = to_string(0) + number;

  WriteXYZ.open("frames/config_" + number + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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