#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);

    ofstream out("es2.csv");
    out << "standard,exponential,lorentz" << endl;

    double lambda = 1.;
    int N = 1e4;
    for (int i = 0; i < N; i++){
        out << rnd.Rannyu() << "," << rnd.Expo(lambda) << "," << rnd.Cauchy(lambda) <<endl;
    }

    out.close();

    return 0;
}
