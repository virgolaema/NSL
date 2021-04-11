#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);

    ofstream out1,out2;
    out1.open("es1.txt");
    out2.open("es2.txt");
    
    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    double * unc = new double [N];
    double * unc_is = new double [N];
    int start = 2; //start from 

    for (int i = start; i < N; i++){ 
              
        int throws = (i+1)*L; //number of throws at this step
        int Nblocks  = i+1; //number of blocks at this step
        cout << "Ho " << Nblocks  << " blocchi, ovvero " << throws  << " valori"<<endl;
        double * A = new double [Nblocks];
        double * A_is = new double [Nblocks];

        for (int j = 0; j < Nblocks ; j++){
            double * appo_A = new double [L];
            double * appo_A_is = new double [L];
            for (int k = 0; k < L; k++){
                appo_A[k] = M_PI/2.*cos(rnd.Rannyu()*M_PI/2.);
                double y_appo = rnd.Rannyu();
                double x_appo =  1 - sqrt(1 - y_appo);
                appo_A_is[k] = M_PI/2. * cos(x_appo*M_PI/2.)*(0.5)*pow(1-x_appo,-1); //normalization 
            }
            A[j] = mean (L, appo_A); 
            A_is[j] = mean (L, appo_A_is); 
            delete []appo_A,appo_A_is;
        }
        
        unc [i] = variance_blocks(Nblocks , A); 
        unc_is [i] = variance_blocks(Nblocks , A_is); 

        out1 << Nblocks << "," << mean (Nblocks,A) << "," << unc[i] << "," << pow(unc[i],0.5) << endl; 
        out2 << Nblocks << "," << mean (Nblocks,A_is) << "," << unc_is[i] << "," << pow(unc_is[i],0.5) << endl; 

        delete []A,A_is;
    }

    out1.close();
    out2.close();
    rnd.SaveSeed();
    return 0;
}
