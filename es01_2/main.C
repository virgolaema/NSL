#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);
    int M = 1e4;                //number of values
    double lambda = 1.;         //parameter for both exp and Cauchy
    int Ni [4] = {1,2,10,100};  //number of elements in the sum

    for (int u = 0; u < 4; u++){
        int N = Ni[u]; //values in each sum SN
        ofstream out("es2_"+to_string(N)+".txt");

        for (int i = 0; i < M; i++){ 
            double appo_un = 0; //uniform distr
            double appo_ex = 0; //expo distr
            double appo_lo = 0; //lorentz distr
            for (int k = 0; k < N; k++){
                appo_un += rnd.Rannyu();
                appo_ex += rnd.Expo(lambda);
                appo_lo += rnd.Cauchy(lambda);
            }
            out << appo_un /(double) N << "," << appo_ex / (double)N << "," << appo_lo  /(double) N << endl; 
        }
        out.close();
    }
    rnd.SaveSeed();
    return 0;
}
