#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    ofstream out1,out2, out3;

    //ex 1.2 and 1.2
    cout << "Exercise 1.1 and 1.2\n";
    out1.open("es1_1.out");
    out2.open("es1_2.out");

    int M = 100000;     //number of gen values
    int N = 100;        // number of blocks
    int L = (int) M/N;  //values per block

    // variables "1" are for the estimation of r, "2" for <(r-0.5)^2>
    double mean1 [N] = {}, mean2 [N] = {};
    double var1 = 0, var2 = 0;
    double sum_var1 = 0., sum_var2 = 0; //progressive sums

    Random rnd;
    rnd.Initialize(rnd);

    for (int i = 0; i < N; i++){ 
        int Nblocks  = i+1;     //number of blocks at this step
        int throws = Nblocks*L; //number of throws at this step
        cout << "Block number " << Nblocks <<endl;

        for (int k = 0; k < L; k++){
            double r = rnd.Rannyu(); 
            mean1[i] += r;
            mean2[i] += pow(r-0.5,2);
        }

        mean1[i] /= (double)L;
        mean2[i] /= (double)L;
        if (Nblocks != 1) { //when N=1, leaves var=0 to avoid -nan 
            var1 = variance_blocks(Nblocks,mean1);
            var2 = variance_blocks(Nblocks,mean2);
        }
        out1 << Nblocks  << "," << mean(Nblocks, mean1) << "," << var1 << endl;
        out2 << Nblocks  << "," << mean(Nblocks, mean2) << "," << var2 << endl;
    }

    out1.close();
    out2.close();

/////////////////////////////////////////////////////////////////////////////////////
    //ex1.1.3
    out3.open("es1_3.out");
    cout << "\n Exercise 1.3\n";

    M = 100;        //number of sub-intervals in 0,1
    int n = 10000;  //gen values in each exp
    int reps = 100; //number of experiments

    double nu = (double) n/M; //expected value per bin
    double chisq [M] = {};
    double nui [M] = {}; //expectation for each bin 
    for (int i = 0; i < M; i++) nui[i] = nu;
    
    for (int k = 0; k < reps; k++){
        cout << "Computing the " << k+1 << " chi-squared\n";
        double ni [M] = {}; //actual values in each bin, initialized to 0
        for (int i = 0; i < n; i++){ //filling the vector ni, repr. the bins
            double r = rnd.Rannyu();
            for (int j = 0; j < M; j++){ 
                if (r >= j*(1./M) and r < (j+1)*1./M) {
                    ni[j]++;
                    j = M; //end cicle, got the correct bin
                }
            }
        }
        out3 << chisquared(M,nui,ni,nui) << endl; //variance is nui, same as exp value
    }

    out3.close();
    rnd.SaveSeed();
    return 0;
}