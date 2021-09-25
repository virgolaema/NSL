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

    //get M (number of generated points)
    int M;
    input >> M;
    //get N (number of blocks)
    int N;
    input >> N;
    //set starting point for wf
    input >> x0;
    //set delta for T and which type of step it is
    input >> delta >> T;
    //params to regulate
    double mu, sigma;
    input >> mu >> sigma;
    //get name for output file
    input >>filename;
    input.close();

    ofstream out1, hist;
    out1.open(filename+".out");
    hist.open("hist.out");

    //create object and set params
    Variational elem (x0);
    elem.setDelta(delta);
    elem.setMu(mu);
    elem.setSigma(sigma);

    double accepted = 0;
    int MRTsteps = 100;     //check if it's okay
    int L = (int) M/N;      //values per block
    int start = 1; //start from the block "start"
    double  A [N] = {};

    for (int i = 0; i < N; i++){ //put N as max
        int Nblocks  = i+1; //number of blocks at this step
        if (Nblocks%5==0) cout <<"--------------------------------------\nComputing with Nblocks = " << Nblocks << endl;
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
                        //accepted++;
                    }
                }
                A[i] += elem.exp_H(elem.X);
                hist << elem.X << "," << pow(elem.phi_T(elem.X),2) << endl;
                //accepted = accepted/(double)MRTsteps*100.; //% rate of acceptance, expected about 50%
                //cout << "Acceptance " << accepted << "%" << endl;
            }
            A[i] /= L;
        }
        out1 << Nblocks << "," << mean(Nblocks,A) << "," << sqrt(variance_blocks(Nblocks, A)) << endl;
    }
    out1.close();
    hist.close();


    //final value of the analysis, need this for grid search for mu, sigma
    ofstream out2;
    out2.open("grids_final.out");
    out2 << mean (N, A);
    out2.close();
    
    //just to display the wavefunction
    ofstream wave("wave.out");
    Variational prova (x0);
    prova.setSigma(sigma);
    prova.setMu(mu);
    double x = -3;
    double step = abs(x)*2/500.;
    for (int i = 0; i < 500; i++){
        x += step;
        prova.setX(x);
        wave << prova.X << "," << pow(prova.phi_T(prova.X),2) << endl;
    }
    wave.close();
    
    elem.rnd.SaveSeed();
    return 0;
}


Variational :: Variational (double x){
    X = x;
    Random appo;
    rnd = appo;
    rnd.Initialize(rnd);
}

Variational :: ~Variational(){}

double Variational:: phi_plus (double x){return exp ( -0.5 * pow((x + mu) / sigma,2));}

double Variational:: phi_minus(double x){return exp ( -0.5 * pow((x - mu) / sigma,2));}

double Variational:: phi_T (double x){return  phi_plus(x) +  phi_minus(x);} 
//normalitazion not relevant in MRT

double Variational:: exp_H (double x) { //<H>
    double kin = -0.5*( phi_minus(x) * pow(x - mu,2)/pow(sigma,4) + phi_plus(x) * pow(x + mu,2)/pow(sigma,4) - phi_T(x)/pow(sigma,2));
    double pot = (pow(x,4) - 2.5 * pow(x,2)) *  phi_T(x);
    return (kin+pot) / phi_T(x) ;
}

double Variational :: A (double xnew, double xold){return min(1., pow(phi_T(xnew) / phi_T(xold),2));}

void Variational :: T (double & xnew, double & xold){xnew = xold + rnd.Rannyu(-delta,delta);}