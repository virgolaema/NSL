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
        double mu = r; 

        void Initialize();
        double d1 (double, double);
        double d2 (double, double);
        double St (double); //discrete case
        double St (double, double); //continuous case
        double N (double);
        double Pt (double, double);
        double Ct (double, double);
};


int main(){

    ofstream out1("Cd.out"), out2("Pd.out"), out3("Cc.out"), out4("Pc.out");
    black_scholes bs;
    bs.Initialize();

    int M = 100000;     //number of gen values
    int N = 100;        //number of blocks
    int L = (int) M/N;  //values per block
    double tfin = bs.T;
    double Ctd [N] = {}, Ptd [N]= {}, Ctc [N]= {}, Ptc [N]= {};

    for (int i = 0; i < N; i++){ 
        int Nblocks  = i+1; //number of blocks at this step
        cout << "Blocks: " << Nblocks << endl;
        double var_Ctd = 0, var_Ptd = 0, var_Ctc = 0, var_Ptc = 0;
        for (int k = 0; k < L; k++){
            //discrete
            Ctd[i] += bs.Ct(bs.St(tfin), tfin);
            Ptd[i] += bs.Pt(bs.St(tfin), tfin);
            
            //continuous
            double t = 0.;
            double Sprec = bs.S0;
            for (int h = 0; h < 100; h++){
                Sprec = bs.St(t, Sprec);
                t += tfin/100.;
            }
            
            Ctc[i] += bs.Ct(bs.St(tfin,Sprec), tfin);
            Ptc[i] += bs.Pt(bs.St(tfin,Sprec), tfin);
        }
        Ctd [i] /= (double) L;
        Ptd [i] /= (double) L;
        Ctc [i] /= (double) L;
        Ptc [i] /= (double) L;
        
        if (Nblocks != 1){ //else, var remains 0
            var_Ctd = variance_blocks(Nblocks, Ctd);
            var_Ptd = variance_blocks(Nblocks, Ptd);
            var_Ctc = variance_blocks(Nblocks, Ctc);
            var_Ptc = variance_blocks(Nblocks, Ptc);
        }
        out1 << Nblocks << "," << mean (Nblocks, Ctd) << "," << var_Ctd << endl;
        out2 << Nblocks << "," << mean (Nblocks, Ptd) << "," << var_Ptd << endl;
        out3 << Nblocks << "," << mean (Nblocks, Ctc) << "," << var_Ctc << endl;
        out4 << Nblocks << "," << mean (Nblocks, Ptc) << "," << var_Ptc << endl;
    }

    out1.close();
    out2.close();
    out3.close();
    out4.close();
    bs.rnd.SaveSeed();
    return 0;
}


void black_scholes :: Initialize (){rnd.Initialize(rnd);}

double black_scholes :: d1 (double S_t, double t){
    return 1./(sigma * sqrt(T-t)) * (log(S_t / K) + (r + pow(sigma,2)/2.) * (T-t));
}

double black_scholes :: d2 (double S_t, double t){
    return d1(S_t,t) - sigma * sqrt(T-t);
}

double black_scholes :: St (double t){
    return S0 * exp ((mu - pow(sigma,2)/2.)*t + sigma * rnd.Gauss(0,t));
}

double black_scholes :: St (double t, double prec){
    return prec * exp (((mu - pow(sigma,2)/2.)*step + sigma * rnd.Gauss(0,1))* sqrt(step));
}

double black_scholes :: N (double x){
    return 0.5 * (1 + erf( x/sqrt(2)));
}

double black_scholes :: Ct (double S_t, double t){
    return S_t * N(d1(S_t,t)) - K * exp(-r*(T-t)) * N(d2(S_t,t));
}

double black_scholes :: Pt (double S_t, double t){
    return S_t * (N(d1(S_t,t)) - 1) - K *exp(-r*(T-t)) * (N(d2(S_t,t)) - 1);
}