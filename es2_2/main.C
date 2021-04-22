#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);

    ofstream out1("es1.txt");
    ofstream out2("es2.txt");

    double a = 1; //dimension of lattice
    double imax = 1e2; //steps till stop
    double x_start = 0, y_start = 0, z_start = 0;
    int M = 1e4; //number of gen values
    int N = 1e1; //values per block
    int L = (int) M/N;
    
    //discrete time
    for (int i = 0; i < imax; i++){ //set steps from 1 to max
        double * A = new double [N];
        for (int j = 0; j < N ; j++){ //move to successive blocks
            double * pos = new double [L];
            for (int k = 0; k < L; k++){ //in the block, make operations
                double x = x_start, y = y_start, z = z_start;
                for (int t = 0; t < i+1; t++){ //steps of the RW

                    double r = rnd.Rannyu();
                    double sign;
                    if (r > 0.5){sign = 1;}
                    else {sign = -1;}

                    double dice = rnd.Rannyu();
                    if (dice >= 1./3. && dice < 2./3.){x+=a*sign;}
                    else if (dice >= 2./3.){y+=a*sign;}
                    else if (dice < 1./3.){z+=a*sign;}
                }        
                pos [k] = x*x + y*y + z*z; //result of single walk in block
            }
            A[j] = mean(L,pos);
            delete []pos;
        }
        double f = sqrt(mean (N,A));
        out1 << i+1 << "," << f << "," << (1./(2*f))*sqrt(variance_blocks(N,A)) << endl; //mean among all blocks and variance
        delete []A;   
    }

    //2nd part
    for (int i = 0; i < imax; i++){ //set steps from 1 to max
        double * A = new double [N];
        for (int j = 0; j < N ; j++){ //move to successive blocks
            double * pos = new double [L];
            for (int k = 0; k < L; k++){ //in the block, make operations
                double x = x_start, y = y_start, z = z_start;
                for (int t = 0; t < i+1; t++){ //steps of the RW

                    double phi = rnd.Rannyu(0,2* M_PI);
                    double theta = rnd.Rannyu(0,M_PI);

                    x += a* cos(phi) * sin(theta);
                    y += a* sin(phi) * sin(theta);
                    z += a* cos(phi);

                }        
                pos [k] = x*x + y*y + z*z; //result of single walk in block
            }
            A[j] = mean(L,pos);
            delete []pos;
        }
        double f = sqrt(mean (N,A));
        out2 << i+1 << "," << f << "," << (1./(2*f))*sqrt(variance_blocks(N,A)) << endl; //mean among all blocks and variance
        delete []A;   
    }

    out2.close();
    out1.close();
    rnd.SaveSeed();
    return 0;
}
