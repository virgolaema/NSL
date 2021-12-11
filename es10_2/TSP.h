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

class Genetic{
    public:
    Genetic(int);
    Random rnd;
    int n_elems;
    int L;
    double p, prob_mut;
    vector <Position> cities; //array of cities
    vector <int> ordered; //needed for check ()
    
    void whichProb (double prob){prob_mut = prob;};
    void whichP (double P){p=P;};
    void whichLoss(int l){L = l;};
    void shuffle (vector <int>&);
    void dispose (int); //0 for circ, 1 for square
    void newIndividual(vector<int> &);
    void check(vector<int>);
    double loss (vector <int>); //module
    void order(vector <Individual> &);
    int selection (vector <Individual>);
    void randomSwap (vector<Individual> &);
    void contSwap (vector<Individual> &);
    void shift (vector<Individual> &);
    void crossover(vector<Individual> &);
    int Pbc (int);
    void end();
};