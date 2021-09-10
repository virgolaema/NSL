#include "../lib.h"

using namespace std;

double hbar = 1, m = 1; //units to make it simpler     

class Variational{ 
    public: 
    double delta; //step of T
    double X, mu, sigma;
    Random rnd;
    Variational (double);
    ~Variational();
    //params
    void setX (double x){X = x;}
    void setDelta (double d){delta = d;} //step of MRT
    void setMu (double m){mu = m;}
    void setSigma(double s){sigma = s;}
    //wavefuncs
    double phi_plus (double);
    double phi_minus(double);
    double phi_T(double);
    double exp_H (double); 
    //MRT
    void T (double&, double&); //Uniform move
    double A (double, double);
};

int main(){
    ifstream input;
    input.open("input.dat");
    double x0, delta;
    int T;
    string filename;

    //set starting point fow wf
    input >> x0;
    //set delta for T and which type of step it is
    input >> delta >> T;
    //get name for output file
    input >>filename;
    input.close();

    ofstream out;
    out.open(filename+".out");

    //create object and set params
    Variational elem (x0);
    elem.setDelta(delta);
    elem.setMu(0.3);
    elem.setSigma(0.1);

    double accepted = 0;
    int MRTsteps = 100; //check if it's okay
    int M = 1000000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    int start = 1; //start from the block "start"

    for (int i = start; i < 2; i++){ //put N as max
        int Nblocks  = i+1; //number of blocks at this step
        cout <<"Computing with Nblocks = " << Nblocks << endl;
        double  A [Nblocks] = {};
        for (int j = 0; j < Nblocks ; j++){
            for (int k = 0; k < L; k++){
                accepted = 0;
                for (int h = 0; h < MRTsteps; h++){
                    double xold = elem.X;
                    double xnew = 0.;
                    //generate x'
                    elem.T(xnew, xold);
                    double alpha = elem.A (xnew, xold);
                    if (elem.rnd.Rannyu() <= alpha) {
                        elem.setX(xnew); //accept the step, update position
                        accepted++;
                    }
                    A[j] += elem.exp_H(elem.X);
                }
                accepted = accepted/(double)MRTsteps*100.; //% rate of acceptance, expected about 50%
            }
            A[j] /= L*(double)MRTsteps;
            cout << "Acceptance " << accepted << "%" << endl;
        }
        out << Nblocks << "," << mean(Nblocks,A) << "," << sqrt(variance_blocks(Nblocks, A)) << endl;
    }
    
    ofstream wave("wave.out");
    Variational prova (0.1);
    prova.setSigma(0.5);
    prova.setMu(0.1);
    double x = -2.2;
    double step = 4.4/500.;
    for (int i = 0; i < 500; i++){
        x += step;
        prova.setX(x);
        wave << prova.X << "," << pow(prova.phi_T(prova.X),2) << endl;
    }
    wave.close();
 

    elem.rnd.SaveSeed();
    out.close();
    return 0;
}


Variational :: Variational (double x){
    X = x;
    Random appo;
    rnd = appo;
    rnd.Initialize(rnd);
}

Variational :: ~Variational(){}

double Variational:: phi_plus (double a){return exp ( -0.5 * pow((a + mu) / sigma,2));}

double Variational:: phi_minus(double a){return exp ( -0.5 * pow((a - mu) / sigma,2));}

double Variational:: phi_T (double a){return  phi_plus(a) +  phi_minus(a);} 
//normalitazion not relevant in MRT

double Variational:: exp_H (double a) { //<H>
    double kin =  ((1 + a) * phi_T(a) +  mu * ( phi_plus(a) -  phi_minus(a)));
    double pot = (pow(a,4) - 2.5 * a) *  phi_T(a);
    return  phi_T(a) * pow (hbar/2./sigma,2) / m * (kin +pot);
}

double Variational :: A (double xnew, double xold){return min(1., pow(phi_T(xnew) / phi_T(xold),2));}

void Variational :: T (double & xnew, double & xold){xnew = xold + rnd.Rannyu(-delta,delta);}