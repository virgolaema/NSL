#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);
    int M = 1e4;
    double lambda = 1.;
    int Ni [4] = {1,2,10,100};

    for (int u = 0; u < 4; u++){
        int N = Ni[u]; //values in each sum SN
        ofstream out("es2_"+to_string(N)+".txt");

        for (int i = 0; i < M; i++){ 

            double * appo_st = new double [N];
            double * appo_ex = new double [N];
            double * appo_lo = new double [N];

            for (int k = 0; k < N; k++){
                appo_st[k] = rnd.Rannyu();
                appo_ex[k] = rnd.Expo(lambda);
                appo_lo[k] = rnd.Cauchy(lambda);
            }
            out << mean (N, appo_st) << "," << mean (N, appo_ex) << "," << mean (N, appo_lo) << endl; 

            delete []appo_st,appo_ex,appo_lo;
        }
            
        out.close();
    }
    rnd.SaveSeed();
    return 0;
}
