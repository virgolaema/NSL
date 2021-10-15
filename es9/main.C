#include "TSP.h"

using namespace std;

int main(){
    int n_cities;  //number of cities
    int n_indiv;  //number of individuals in each gen
    int n_gens;   //number of generations
    int disp, loss, p;

    ifstream in ("input.dat");
    in >> n_cities;
    cout << "Number of cities " << n_cities << endl;
    in >> n_indiv;
    in >> n_gens;
    cout << "Number of generations " << n_gens << endl;
    in >> disp;
    in >> loss;
    cout << "Loss function " << loss << endl;
    in >> p;
    in.close();

    Genetic tsp (n_cities);
    tsp.dispose(disp);        //dispose the n cities on a circ (0) or a square (1)
    tsp.whichLoss(loss);      //1 for L1, 2 for L2
    tsp.whichP(p);         //exp of selection operator
    tsp.whichProb(0.1);    //probabilty of mutations

    ofstream pos ("cities.out");
    for (int i = 0; i < n_cities; i++) pos << tsp.cities[i].getX() << "," << tsp.cities[i].getY() << endl;
    pos.close();

    //Create first generation
    vector <Individual> generation; 
    generation.resize(n_indiv);
    for (int i = 0; i < n_indiv; i++){
        generation[i].setN(n_cities);
        tsp.newIndividual(generation[i].sequence); //start creating random individuals
        tsp.check(generation[i].sequence);
        generation[i].fitness = tsp.loss(generation[i].sequence);
    }
    tsp.order(generation);

    //Create new generations
    ofstream fit ("fitness.out");
    ofstream meanL ("meanL.out");

    for (int i = 1; i <= n_gens; i++){ //should be n_gens
        if(i%10 == 1) cout << "Creating the " << i << " generation" << endl;
        //NEW GEN
        tsp.crossover(generation);
        //MUTATIONS. Every mutation also puts element in fitness order
        tsp.randomSwap(generation); 
        tsp.contSwap(generation);
        tsp.shift(generation);
        if(i%10 ==0) cout << "Loss of best element at this gen: " << generation[0].fitness << endl;
        fit << generation[0].fitness << endl;
        double mean_L = 0.;
        for (int h = 0; h < n_indiv/2.; h++){
            mean_L += generation[h].fitness;
        }
        meanL << mean_L/(n_indiv/2.) << endl;
    }
    meanL.close();
    fit.close();

    ofstream final("final.out");
    cout << "Printing best configuration...\n"; 
    for (int h = 0; h < n_cities; h++) final << tsp.cities[generation[0].sequence[h]].getX() << "," << tsp.cities[generation[0].sequence[h]].getY() << endl;
    final << tsp.cities[generation[0].sequence[0]].getX() << "," << tsp.cities[generation[0].sequence[0]].getY() << endl;
    final.close();

    return 0;
}