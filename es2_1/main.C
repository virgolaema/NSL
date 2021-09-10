#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){
    Random rnd;
    rnd.Initialize(rnd);

    ofstream out1,out2;
    out1.open("es1.out");
    out2.open("es2.out");
    
    int M = 100000;         //number of gen values
    int N = 100;            //number of blocks
    int L = (int) M/N;      //values per block
    double A_un [N] = {}; //uniform
    double A_is [N] = {}; //importance sampling
    double var_un = 0, var_is = 0;
    double x_appo, y_appo;

    for (int i = 0; i < N; i++){ 
    
        int Nblocks  = i+1;     //number of blocks at this step
        int throws = Nblocks*L;   //number of throws at this step
        cout << "Block " << Nblocks <<endl;
        
        for (int k = 0; k < L; k++){
            // sampling with uniform probability
            A_un[i] += M_PI/2.*cos(rnd.Rannyu()*M_PI/2.);
            // importance sampling
            y_appo = rnd.Rannyu();
            x_appo = 1 - sqrt(1 - y_appo);
            A_is[i] += M_PI/2. * cos(x_appo*M_PI/2.)/(2-2*x_appo); //problem with constants
        }

        A_un[i] /= (double) L; 
        A_is[i] /= (double) L; 
        
        if (Nblocks != 1){ //avoid -nan when Nblocks = 1, var remains 0 as initialization
            var_un = variance_blocks(Nblocks , A_un); 
            var_is = variance_blocks(Nblocks , A_is); 
        }

        out1 << Nblocks << "," << mean(Nblocks, A_un) << "," << var_un << endl; 
        out2 << Nblocks << "," << mean(Nblocks, A_is) << "," << var_is << endl; 
    }

    out1.close();
    out2.close();
    rnd.SaveSeed();

    return 0;
}
