#include "../lib.h"
#include "../random.h"

using namespace std;

class black_scholes{
    public:
        Random rnd;
        double S0 = 100;
        double T = 1.;
        double K = 100.;
        double r = 0.1;
        double sigma = 0.25;
        double step = 0.01; //step in continuous case
        double mu = r; //drift, correct?

        void Initialize();
        double d1 (double, double);
        double d2 (double, double);
        double St (double); //discrete case
        double St (double, double); //continuous case
        double N (double);
        double Pt (double, double);
        double Ct (double, double);

    private:
};

void black_scholes :: Initialize (){rnd.Initialize(rnd);}
double black_scholes :: d1 (double S_t, double t){return 1./(sigma * sqrt(T-t)) * (log(S_t / K) + (r + pow(sigma,2)/2.) * (T-t));}
double black_scholes :: d2 (double S_t, double t){return d1(S_t,t) - sigma * sqrt(T-t);}
double black_scholes :: St (double t){return S0 * exp ((mu - pow(sigma,2)/2.)*t + sigma * rnd.Gauss(0,t));}
double black_scholes :: St (double t, double prec){return prec * exp (((mu - pow(sigma,2)/2.)*step + sigma * rnd.Gauss(0,1))* sqrt(step));}
double black_scholes :: N (double x){return 0.5 * (1 + erf( x/sqrt(2)));}
double black_scholes :: Ct (double S_t, double t){return S_t * N(d1(S_t,t)) - K * exp(-r*(T-t)) * N(d2(S_t,t));}
double black_scholes :: Pt (double S_t, double t){return S_t * (N(d1(S_t,t)) - 1) - K *exp(-r*(T-t)) * (N(d2(S_t,t)) - 1);}

int main(){

    ofstream out1("Cd.txt"), out2("Pd.txt"), out3("Cc.txt"), out4("Pc.txt");
    black_scholes bs;
    bs.Initialize();

    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    int start = 1; //start from the block "start"
    double tfin = bs.T;
//N=3;
    for (int i = start; i < N; i++){ 
        int Nblocks  = i+1; //number of blocks at this step
        double Ctd [Nblocks] = {}, Ptd [Nblocks]= {}, Ctc [Nblocks]= {}, Ptc [Nblocks]= {};
        for (int j = 0; j < Nblocks ; j++){
            for (int k = 0; k < L; k++){
                //discrete
                Ctd [j] += bs.Ct(bs.St(tfin), tfin);
                Ptd [j] += bs.Pt(bs.St(tfin), tfin);
                //cout << "S finale discreto " << bs.St(tfin) << endl;
                //continuous
                double t = 0.;
                double Sprec = bs.S0;
                for (int h = 0; h < 100; h++){
                    Sprec = bs.St(t, Sprec);
                    t += tfin/100.;
                    //cout << "Tempo " << t << endl;
                }
                //cout << "S finale continuo " << Sprec << endl;
                
                Ctc [j] += bs.Ct(bs.St(tfin,Sprec), tfin);
                Ptc [j] += bs.Pt(bs.St(tfin,Sprec), tfin);
            }
            Ctd [j] /= L;
            Ptd [j] /= L;
            Ctc [j] /= L;
            Ptc [j] /= L;
        }
        out1 << Nblocks << "," << mean (Nblocks, Ctd) << "," << sqrt(variance_blocks(Nblocks, Ctd)) << endl;
        out2 << Nblocks << "," << mean (Nblocks, Ptd) << "," << sqrt(variance_blocks(Nblocks, Ptd)) << endl;
        out3 << Nblocks << "," << mean (Nblocks, Ctc) << "," << sqrt(variance_blocks(Nblocks, Ctc)) << endl;
        out4 << Nblocks << "," << mean (Nblocks, Ptc) << "," << sqrt(variance_blocks(Nblocks, Ptc)) << endl;
        cout << "step " << i << endl;

    }

    out1.close();
    out2.close();
    out3.close();
    out4.close();
    bs.rnd.SaveSeed();
    return 0;
}
