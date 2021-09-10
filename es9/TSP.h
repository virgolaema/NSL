#include "../lib.h"
#include "../random.h"

using namespace std;

class Individual{
    public:
    int N;
    vector <int> sequence;
    double fitness;
    void setN(int); //instead of a initializer
};

class Genetic{
    public:
    Genetic(int);
    Random rnd;
    int n_elems;
    int L;
    Position * cities; //array of cities
    vector <int> ordered; //needed for check ()

    void whichLoss(int l){L = l;};
    void shuffle (vector <int>&);
    void dispose (int); //0 for circ, 1 for square
    void newIndividual(vector<int> &);
    void check(vector<int>);
    double loss (vector <int>); //module
    void order(vector <Individual> &);
    void selection (vector <Individual>, Individual&);
    void newGen (vector <Individual>, vector <Individual> &);
    void randomSwap (vector<Individual> &);
    void shift (vector<Individual> &);
    void end();
};