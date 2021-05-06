#include "../lib.h"
#include "../random.h"

using namespace std;

class Hydrogen: public Position{
    public: 
    double delta; //step of T
    double r_expected;
    int qns; //quantum numbers, ie 100 or 210
    Random rnd;

    Hydrogen (double, double, double);
    ~Hydrogen();
    void setQns(int);
    void setDelta (double d){delta = d;}
    double getDelta(){return delta;}
    double getRexp (){return r_expected;}
    double p (Position &);
    void TU (Position&, Position&); //Uniform move
    void TG (Position &pnew, Position &pold); //Gaussian move
    double q (Position &, Position &);
    double A (Position&, Position&);
    protected:
};

int main(){
    ifstream input;
    input.open("input.dat");
    double x0, y0, z0, delta;
    int qns, T;
    string filename;

    //set starting point fow wf
    input >> x0 >> y0 >> z0;
    //set quantum numbers 
    input >> qns;
    //set delta for T and which type of step it is
    input >> delta >> T;
    //get name for output file
    input >>filename;
    input.close();

    ofstream out;
    out.open(filename+".out");

    //create object and set params
    Hydrogen H (x0, y0, z0);
    H.setQns(qns);
    H.setDelta(delta);

    double accepted = 0;
    int MRTsteps = 100; //check if it's okay
    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    int start = 2; //start from the block "start"

    for (int i = start; i < N; i++){ 
        int Nblocks  = i+1; //number of blocks at this step
        cout <<"Computing with Nblocks = " << Nblocks << endl;
        double  A [Nblocks] = {};
        for (int j = 0; j < Nblocks ; j++){
            accepted = 0;
            for (int k = 0; k < L; k++){
                for (int h = 0; h < MRTsteps; h++){
                    Position pold (H.getX(), H.getY(), H.getZ());
                    Position pnew(0.,0.,0.);
                    //generate x'
                    if(T == 0){H.TU(pnew, pold);} 
                    else if (T == 1){H.TG(pnew,pold);}

                    double alpha = H.A (pnew, pold);
                    if (H.rnd.Rannyu() <= alpha) {
                        H.copyPos(pnew); //accept the step, update position
                        accepted++;
                    }
                }
                accepted /= 100./MRTsteps; //% rate of acceptance, expected about 50%
                A[j] += H.R();
            }
            A[j] /= L;
        }
        out << Nblocks << "," << mean(Nblocks,A) - H.getRexp() << "," << sqrt(variance_blocks(Nblocks, A)) << endl;
    }

    H.rnd.SaveSeed();
    out.close();
    return 0;
}


Hydrogen :: Hydrogen (double x, double y, double z):Position (x,y,z){
    Random appo;
    rnd = appo;
    rnd.Initialize(rnd);
}

Hydrogen :: ~Hydrogen(){}

void Hydrogen :: setQns (int q){
    qns = q; //set quantum numbers, which wavefunc between 100 and 210
    if(qns == 100){r_expected = 3./2.;}  //a0 units
    else if (qns == 210){r_expected = 5.;}
}

double Hydrogen :: p (Position& a){ 
    double prob = 0;
    if (qns == 210) {prob = a.getX()*a.getX()/M_PI/32. * exp(-a.R());}
    else if (qns == 100) {prob = exp(-2.*a.R()) / M_PI;}
    return prob;
}

double Hydrogen :: q (Position& pnew, Position& pold){return p(pnew) / p(pold);} 

double Hydrogen :: A (Position& pnew, Position& pold){return min(1., q(pnew,pold));}

void Hydrogen :: TU (Position &pnew, Position &pold){
    pnew.setX(pold.getX() + rnd.Rannyu(-delta,delta));
    pnew.setY(pold.getY() + rnd.Rannyu(-delta,delta)); 
    pnew.setZ(pold.getZ() + rnd.Rannyu(-delta,delta)); 
}

void Hydrogen :: TG (Position &pnew, Position &pold){
    pnew.setX(pold.getX() + rnd.Gauss(0,delta));
    pnew.setY(pold.getY() + rnd.Gauss(0,delta)); 
    pnew.setZ(pold.getZ() + rnd.Gauss(0,delta)); 
}
