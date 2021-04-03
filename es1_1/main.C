#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    ofstream out;
    out.open("es1_1.csv");
    out << "Nblocks,throws,mean_r,var_r,stdev_r,mean_var,stdev_var" << endl;

    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    double * unc = new double [N];
    double * r = new double [N];
    int start = 9; //start from 4 blocks

    Random rnd;
    rnd.Initialize(rnd);

    for (int i = start; i < N; i++){ 
                
        int throws = (i+1)*L; //number of throws at this step
        int Nblocks  = i+1; //number of blocks at this step
        cout << "Ho " << Nblocks  << " blocchi, ovvero " << throws  << " valori"<<endl;
        double * A = new double [Nblocks];
        double * varA = new double [Nblocks];

        for (int j = 0; j < Nblocks ; j++){
            double * appo_A = new double [L];
            double * appo_varA = new double [L];

            //cout << "creato vettore ausiliario di " << L << " elementi "<<endl;
            for (int k = 0; k < L; k++){
                appo_A[k] = rnd.Rannyu();
                appo_varA[k] = pow(appo_A[k] - 0.5,2);  
            }
            A[j] = mean (L, appo_A); 
            varA[j] = mean (L, appo_varA); 

            delete []appo_A,appo_varA;
        }
        
        unc [i] = variance_blocks(Nblocks , A); 

        //"Nblocks,throws,mean_r,var_r,stdev_r,mean_var,stdev_var"
        out << Nblocks  << "," << throws << "," << mean (Nblocks,A) << "," << unc[i] << "," << pow(unc[i],0.5) 
            << "," << mean(Nblocks,varA) << "," << pow(variance_blocks(Nblocks ,varA),0.5) << endl; //save all in csv file 
        
        delete []A,varA;
        //in.close();
    }

    delete []r,unc;
    out.close();

/////////////////////////////////////////////////////////////////////////////////////////////////////7
    //ex1.1.3

    M = 100; //number of sub-intervals in 0,1
    int n = 10000;
    int reps = 100;

    double nu = (double) n/M; //expected value per bin
    double * chisq = new double [M];
    double * nui = new double [M]; //values expected in every bin (and also variance), need this for chisquared function
    for (int i = 0; i < M; i++){nui[i] = nu;}

    out.open("es1_3.csv");

    out << "chi_squared" << endl;
    
    for (int k = 0; k < reps; k++){
        cout << "Calcolo il " << k+1 << " chi-quadro " << endl;
        double * ni = new double [M];
        for (int i = 0; i < M; i++){ni[i] = 0;} //initialize to 0 all the elements

        double appo = 0;
        for (int i = 0; i < n; i++){ //filling the vector ni, repr. the bins
            in >> appo;
            for (int j = 0; j < M; j++){
                if (appo >= j*(1./M) and appo < (j+1)*1./M) {
                    ni[j]++;
                    j = M; //end cicle, got the correct bin
                }
            }
        }
 
        cout << "Chi quadro vale " << chisquared(M,nui,ni,nui) << endl;
        out << chisquared(M,nui,ni,nui) << endl;
        delete []ni;
    }

    delete []chisq;
    out.close();
    in.close();

    rnd.SaveSeed();
    return 0;
}
