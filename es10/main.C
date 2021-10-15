#include "TSP.h"

using namespace std;

int main(){
    int n_cities;  //number of cities
    int blk_steps, MC_Steps;    //max number of steps
    int disp, loss;

    ifstream in ("input.dat");
    in >> n_cities;
    cout << "Number of cities " << n_cities << endl;
    in >> blk_steps;
    in >> MC_Steps;
    in >> disp;
    in >> loss;
    cout << "Loss function " << loss << endl;
    in.close();

    SimAnn tsp (n_cities);
    tsp.dispose(disp);        //dispose the n cities on a circ (0) or a square (1)
    tsp.whichLoss(loss);      //1 for L1, 2 for L2
    tsp.whichProb(1);    //probabilty of mutations

    ofstream pos ("cities.out");
    for (int i = 0; i < n_cities; i++) pos << tsp.cities[i].getX() << "," << tsp.cities[i].getY() << endl;
    pos.close();

    Individual old, attempt;
    old.setN(n_cities);
    tsp.newIndividual(old.sequence); //start creating random individuals
    old.fitness = tsp.loss(old.sequence);

    double Tin = 1000., T = Tin, P;
    double beta = 1./T;
    int count = 0;
    ofstream fit ("fitness.out");
    for(int i=0; i< MC_Steps; i++){
        cout << "MC step " << i << endl;
        cout << "beta " << beta << " T " << T << endl;
        for (int h = 0; h < blk_steps; h++){
            attempt = old;
            //MUTATIONS (compute fitness after every mutation)
            int which = tsp.rnd.Rannyu(1,4);
            if (which == 1) tsp.randomSwap(attempt); 
            if (which == 2) tsp.contSwap(attempt);
            if (which == 3) tsp.shift(attempt);
            P = exp(-(beta)*(attempt.fitness-old.fitness));
            //else if (attempt.fitness < old.fitness){P = 1.;};
            if (tsp.rnd.Rannyu() < P) {
                //cout << "fitness old " << old.fitness << " fitness new " << attempt.fitness << endl;
                old = attempt;
                count++;
            }
        }
        fit << old.fitness << endl;
        cout << old.fitness << endl;
        T *= 0.8; //divide by 2 at each step
        beta = 1./T;
//        if (count < 5) break; //if the element doesn't change anymore, we're at the gs
cout << "numero scambi " << count << endl;
        count = 0;
    }
    fit.close();

    ofstream final("final.out");
    cout << "Printing best configuration...\n"; 
    for (int h = 0; h < n_cities; h++) final << tsp.cities[old.sequence[h]].getX() << "," << tsp.cities[old.sequence[h]].getY() << endl;
    final << tsp.cities[old.sequence[0]].getX() << "," << tsp.cities[old.sequence[0]].getY() << endl;
    final.close();

    return 0;
}