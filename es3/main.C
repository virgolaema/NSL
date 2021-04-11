#include "../lib.h"
#include "../random.h"

using namespace std;

Random rnd;

//defining these as global to ease the function definition
double S0 = 100;
double T = 1;
double K = 100;
double r = 0.1;
double sigma = 0.25;
double mu = 0; //drift, what is it?


double d1 (double);
double d1 (double t){
    return 1./(sigma * sqrt(T-t)) * (log(S0 / K) + (r + pow(sigma,2)/2.) * (T-t));
}

double d2 (double);
double d2 (double t){
    return d1(t) - sigma * sqrt(T-t);
}

double St (double);
double St (double t){
    return S0 * exp ((mu - pow(sigma,2)/2.)*t + sigma * rnd.Gauss(0,t));
}

double N (double);
double N (double x){
    return 0.5*(1+erf(x/sqrt(2)));
}

double Ct (double);
double Ct (double t){ 
    return St(t) * N(d1(t)) - K * exp(-r*(T-t)) * N(d2(t));
}

double Pt (double);
double Pt (double t){ 
    return St(t) * (N(d1(t)) - 1) - K *exp(-r*(T-t)) * (N(d2(t)) - 1);
}

int main(){

    ofstream out1("es1.txt");
    rnd.Initialize(rnd);

    double t0 = 0;

    double t=t0;
    cout << "t=0 " << Ct(t) << "," << Pt(t) << endl;
    cout << "t=T " << Ct(T) << "," << Pt(T) << endl;


    out1.close();
    rnd.SaveSeed();
    return 0;
}
