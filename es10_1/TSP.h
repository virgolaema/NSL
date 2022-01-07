#include "../lib.h"
#include "../random.h"

using namespace std;

class Individual{
    public:
    int N;
    vector <int> sequence;
    double fitness;
    void setN(int n){N = n;}; //instead of a initializer
    private:
};

class SimAnn{
    public:
    SimAnn(int);
    Random rnd;
    int n_elems;
    int L;
    double p, prob_mut, dbeta = 1e-6; //check this
    vector <Position> cities; //array of cities
    vector <int> ordered; //needed for check ()
    
    void whichProb (double prob){prob_mut = prob;};
    void whichLoss(int l){L = l;};
    void shuffle (vector <int>&);
    void dispose (int); //0 for circ, 1 for square
    void newIndividual(vector<int> &);
    void check(vector<int>);
    double loss (vector <int>); //module
    void order(vector <Individual> &);
    void randomSwap (Individual &);
    void contSwap (Individual &);
    void shift (Individual &);
    int Pbc (int);
    void end();
};