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
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  for(int iblk=1;iblk<=nblk;iblk++){
     Reset(iblk);
     cout << "Block " << iblk << endl;
      for(int istep=0;istep<nstep/10;istep++){
        Move();
        Measure();
        Accumulate();
      }
    Average(iblk);
  }
  ConfFinal();  
  return 0;  
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;
  int oldconf, rescale;

  // cout << "Do you want to use the configuration at step -dt?" << endl;
  // cout << "Type 1 to confirm, 0 to decline:" << endl;  
  // cin >> oldconf;
  oldconf = 1;
  // cout << "Do you want to rescale the coordinates in order to match the desired T?" << endl;
  // cout << "Type 1 to confirm, 0 to decline:" << endl;
  // cin >> rescale;
  rescale = 1;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double) npart/rho;
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
  iw = 1;
  n_props = 4; //Number of observables
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){
    vx[i] = rand()/double(RAND_MAX) - 0.5;
    vy[i] = rand()/double(RAND_MAX) - 0.5;
    vz[i] = rand()/double(RAND_MAX) - 0.5;

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
  double sumv2 = 0.0, fs;
  for (int i=0; i<npart; ++i){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;

  if (oldconf==0) {
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  else{
    ReadConf.open("config.00");
    for (int i=0; i<npart; ++i){
    ReadConf >> xold[i];
    ReadConf >> yold[i];
    ReadConf >> zold[i];
   }
    
  }

  ReadConf.close();

  if (rescale == 1){
    double v2tot = 0;
    for (int i=0; i<npart; ++i){
      vx[i] = (x[i] - xold[i])/delta;
      vy[i] = (y[i] - yold[i])/delta;
      vz[i] = (z[i] - zold[i])/delta;

      v2tot += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]; 
    }
    double v2tot_new = 3*temp*npart;
    double cf = v2tot_new / v2tot; 
    for (int i=0; i<npart; ++i){
      xold[i] = Pbc(x[i] - sqrt(cf)* vx[i] *delta);
      yold[i] = Pbc(y[i] - sqrt(cf)* vy[i] *delta);
      zold[i] = Pbc(z[i] - sqrt(cf)* vz[i] *delta);
    }
  } 
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

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
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
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij,wij,w;
  double dx, dy, dz, dr;
  ofstream Epot, Pres;
  v = 0.0; //Reset observables
  t = 0.0;
  w= 0.0;
  //cycle over pairs of particles
  for (int k=0; k<nbins; ++k)  hist_g[k]=0;
    for (int i=0; i<npart-1; ++i){
      for (int j=i+1; j<npart; ++j){

      dx = Pbc( xold[i] - xold[j] ); 
      dy = Pbc( yold[i] - yold[j] ); 
      dz = Pbc( zold[i] - zold[j] ); 

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      //update histogram of g(r)
  
      for (int k=0; k<nbins; ++k){
        if(dr>bin_size*k and dr<=bin_size*(k+1)){
            hist_g[k]+=2;
            break;
        }
      }

      if(dr < rcut){
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
        //Potential energy
        v += vij;
        w += wij;
        
      }          
    }
  }
  //Kinetic energy
  
  walker_av[iv] = v;  // /(double)npart; //Potential energy per particle
  walker_av[iw] = (48.0 * w / 3.0);// /(double)npart;

  // Epot << walker_av[iv]  << endl;
  // Pres << walker_av[iw] << endl;
  return;
}

void Accumulate(){
  for(int i=0; i<nbins;i++){
    hist_ave[i]=hist_ave[i]+hist_g[i];
  }
  for(int i=0; i<2; i++){
    blk_av[i]=blk_av[i]+walker_av[i];
  }
  norm=norm+1;
}

void Reset(int iblk){
  if( iblk==1){
      for(unsigned int i =0; i<nbins; i++){
        glob_av2[i]=0;
        glob_av[i]=0;
      }
    }
  for(unsigned int i=0;i < nbins; i++){
    hist_ave[i]=0;
  }
  for(unsigned int i=0; i<2 ; i++){
    blk_av[i]=0;
  }
  norm=0;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;

  WriteConf.open("config.final");
  for (int i=0; i<npart; ++i)    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  WriteConf.close();

  WriteConf.open("old.final");
  for (int i=0; i<npart; ++i) WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  WriteConf.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
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

void block_error(double* error, double* sum_prog, double* x1, double* x2, int n){
   double* sum_prog2= new double[n];
   
   for (int i=0; i<n; i++){
      sum_prog[i]=0;
		for (int j=0; j<i+1; j++){
         
			sum_prog[i] += x1[j];
         sum_prog2[i] +=x2[j];

		}
		
		sum_prog[i] /= (i+1);
      sum_prog2[i] /= (i+1);
	}

   for(int i=0; i<n; i++){
      error[i]= sqrt( (sum_prog2[i]-pow(sum_prog[i],2)) / i);

   }
  delete[] sum_prog2;

}

void BlockStat(string filein, string fileout){

    int N=100;
    int M=1000;
    double* x1=new double[N]();
    double* x2=new double[N]();
    double* sum_prog=new double[N]();
    double* error=new double[N]();
    int L=M/N;
    ifstream in;
    ofstream out;
    in.open(filein);

    for(int i=0; i < N; i++){
        double appo1=0;
        double sum=0;
        for(int k=0; k < L; k++){
            in >> appo1;
            sum+=appo1;
        }
        x1[i]=sum/L;
        x2[i]=pow(x1[i],2);
    }

    block_error(error,sum_prog,x1,x2,N);
    in.close();

    out.open(fileout);
    for(int i=0; i < N; i++){

      out << sum_prog[i] << " " << error[i] << endl;

    }

  out.close();
  delete[] x1;
  delete[] x2;
  delete[] sum_prog;
  delete[] error;
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Average(int iblk){
  double err_pot, err_press;
  ofstream Gave,epot,pres;
  epot.open("ave_epot.out",ios::app);
  pres.open("ave_pres.out",ios::app);
  Gave.open("output_gave.out",ios::app);
  double r=0;
  double DV=0;
  const int wd=12;

    stima_pot = blk_av[iv]/norm/(double)npart;  //Potential energy
    global_m[iv] += stima_pot;
    global_m2[iv] += stima_pot*stima_pot;
    err_pot=Error(global_m[iv],global_m2[iv],iblk);

    stima_pres = rho * temp + (blk_av[iw]/norm  +  (double)npart) / vol; //Pressure
    global_m[iw] += stima_pres;
    global_m2[iw] += stima_pres*stima_pres;
    err_press=Error(global_m[iw],global_m2[iw],iblk);
    
    //Potential energy per particle
    epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << global_m[iv]/(double)iblk << setw(wd) << err_pot << endl;
    //Pressure
    pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << global_m[iw]/(double)iblk << setw(wd) << err_press << endl;
    
    for(unsigned int i=0; i<nbins; i++){
        r=i*bin_size;
        DV = (4./3.)*M_PI*(pow(r+bin_size,3)-pow(r,3));
        stima_g = (double)hist_ave[i]/(npart*DV*rho*norm);
        glob_av[i] += stima_g;
        glob_av2[i] += stima_g*stima_g;
        err_gdir[i]=Error(glob_av[i], glob_av2[i],iblk);
    }

    if( iblk==nblk){
      for(unsigned int i=0; i<nbins; i++){
          Gave << glob_av[i]/((double)iblk) << " " << err_gdir[i] << endl;
      }
    }
    Gave.close();
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
