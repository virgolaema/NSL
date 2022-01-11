#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){

    Random rnd;
    rnd.Initialize(rnd);

    ofstream out1("es1.txt");
    ofstream out2("es2.txt");

    double a = 1;       //dimension of lattice
    double Nmax = 100;  //max number of steps
    double x_start = 0, y_start = 0, z_start = 0;
    double x, y, z;
    int M = 10000;      //number of walks
    
    //discrete steps
    cout << "Discrete random walk" << endl;
    for (int N = 1; N <= Nmax; N++){ //set time steps from 1 to Nmax
        cout << "Time steps: " << N << endl;
        double A [M] = {};
        for (int j = 0; j < M; j++){ //move to successive blocks
            
            x = x_start; 
            y = y_start; 
            z = z_start;
            
            for (int t = 0; t < N+1; t++){ //time steps of the RW
                double r = rnd.Rannyu();
                double sign;
                if (r > 0.5){sign = 1;}
                else {sign = -1;}

                double dice = rnd.Rannyu();
                if (dice >= 1./3. && dice < 2./3.){x += a*sign;}
                else if (dice >= 2./3.){y += a*sign;}
                else if (dice < 1./3.){z += a*sign;}
                    
                A[j] = x*x + y*y + z*z; //result of single walk 
            }    
        }
        double mean_pos = sqrt(mean (M,A));
        out1 << N << "," << mean_pos << "," << (1./(2*mean_pos))*sqrt(variance(M,A)) << endl; 
    }

    //continuum
    cout << "\nContinuous random walk\n";
    for (int N = 1; N <= Nmax; N++){ //set time steps from 1 to max
        cout << "Time steps: " << N << endl;
        double A [M] = {};
        for (int j = 0; j < M ; j++){ 
            x = x_start; 
            y = y_start; 
            z = z_start;
            
            for (int t = 0; t < N+1; t++){ //steps of the RW

                double phi = rnd.Rannyu(0,2* M_PI);
                double theta = rnd.Rannyu(0,M_PI);

                x += a* cos(phi) * sin(theta);
                y += a* sin(phi) * sin(theta);
                z += a* cos(phi);

            }        
            A[j] = x*x + y*y + z*z; //result of single walk in block
        }
        double mean_pos = sqrt(mean (M,A));
        out2 << N << "," << mean_pos << "," << (1./(2*mean_pos))*sqrt(variance(M,A)) << endl;   
    }

    out2.close();
    out1.close();
    rnd.SaveSeed();
    return 0;
}
