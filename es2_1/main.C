#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);

    ofstream out;
    out.open("es2_1.csv");
    out << "Nblocks,throws,mean_I,var_I,stdev_I" << endl;

    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    double * unc = new double [N];
    double * r = new double [N];
    int start = 2; //start from 

    for (int i = start; i < N; i++){ 
              
        int throws = (i+1)*L; //number of throws at this step
        int Nblocks  = i+1; //number of blocks at this step
        cout << "Ho " << Nblocks  << " blocchi, ovvero " << throws  << " valori"<<endl;
        double * A = new double [Nblocks];

        for (int j = 0; j < Nblocks ; j++){
            double * appo_A = new double [L];
            for (int k = 0; k < L; k++){
                appo_A[k] = cos(rnd.Rannyu()*M_PI/2.);
            }
            A[j] = mean (L, appo_A) * M_PI/2; 

            delete []appo_A;
        }
        
        unc [i] = variance_blocks(Nblocks , A); 

        //"Nblocks,throws,I,var_I,stdev_I"
        out << Nblocks  << "," << throws << "," << mean (Nblocks,A) << "," << unc[i] << "," << pow(unc[i],0.5) << endl; //save all in csv file 
        
        delete []A;
        //in.close();
    }



    rnd.SaveSeed();
    return 0;
}
