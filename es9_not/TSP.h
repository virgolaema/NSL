#include "../lib.h"
#include "../random.h"

using namespace std;

class Genetic{
    public:
    Ising2D(int);
    Random rnd;
    int n_elems;
    int ** matrix;
    Position * cities; 

    void set1 (int, int);
    void set0 (int, int);
    void dispose (int); //0 for circ, 1 for square
    void swap (int, int,int, int);
    void check();
    double loss1 (); //module
    double loss2 (); //squared
    void end();
};

