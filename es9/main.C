#include "TSP.h"

using namespace std;

int main(){
    int n_cities = 32; //number of cities
    int n_indiv = 500; //number of individuals in each gen
    int n_gens = 500; //number of generations

    Genetic circ (n_cities);
    circ.dispose(0);        //dispose the n cities on a circumference
    circ.whichLoss(1);      //1 for L1, 2 for L2
    circ.whichP(2);         //exp of selection operator
    circ.whichProb(0.8);    //probabilty of mutations

    ofstream pos ("circ.out");
    for (int i = 0; i < n_cities; i++) pos << circ.cities[i].getX() << "," << circ.cities[i].getY() << endl;
    pos.close();
    cout << "Elements disposed on a circ \n";

    vector <Individual> previous_gen; 
    previous_gen.resize(n_indiv);
    for (int i = 0; i < n_indiv; i++){
        previous_gen[i].setN(n_cities);
        circ.newIndividual(previous_gen[i].sequence); //start creating random individuals
        circ.check(previous_gen[i].sequence);
        previous_gen[i].fitness = circ.loss(previous_gen[i].sequence);
    }
    circ.order(previous_gen);

    // creating a vector of individuals, that will be updated with the next gen everytime
    vector <Individual> next_gen = previous_gen;

    //NEW GENS
    ofstream fit ("fitness.out");
    for (int i = 1; i < n_gens; i++){ //should be n_gens
        cout << "Creating the " << i+1 << " generation..."<< endl;
        //MUTATIONS. Every mutation also puts element in fitness order
        circ.crossover(previous_gen);
        circ.randomSwap(previous_gen); 
        circ.contSwap(previous_gen);
        circ.shift(previous_gen);
        //mutations done, now create next gen
        circ.newGen(previous_gen, next_gen);
        circ.order(next_gen);
        previous_gen = next_gen;
        cout << "loss of first element, should be the best: " << next_gen[0].fitness << endl;
        fit << next_gen[0].fitness << endl;
    }
    fit.close();

    for (int i = 0; i < n_indiv; i++){
        //cout << "Fitness of " << i << " element is " << next_gen[i].fitness << endl;
    }


    /*for (int i = 0; i < n; i++){
        cout << "\n The " << i+1 <<" element of the last created generation is: \n";
        for (int h = 0; h < n; h++) cout << next_gen[i].sequence[h] << " ";
        cout << endl;
    }*/

    ofstream final("final.out");

    cout << "printing best configuration..\n"; 
    for (int h = 0; h < n_cities; h++) 
        final << circ.cities[next_gen[0].sequence[h]-1].getX() << "," << circ.cities[next_gen[0].sequence[h]-1].getY() << endl;
    final << circ.cities[next_gen[0].sequence[0]-1].getX() << "," << circ.cities[next_gen[0].sequence[0]-1].getY() << endl;

    final.close();

    //cout << "loss of first element, should be the best: " << next_gen[0].fitness << endl;

    //now, get busy. How many generations? selection mutation

    return 0;
}

 /*ofstream disp_circ ("circ.out");
    ofstream disp_sqre ("sqre.out");

    for (int i = 0; i < n; i++){
        disp_circ << circ.cities[i].getX() << "," << circ.cities[i].getY() << endl;
        disp_sqre << sqre.cities[i].getX() << "," << sqre.cities[i].getY() << endl;
    }

    disp_sqre.close();
    disp_circ.close();

    sqre.end();
    circ.end();*/

    /*
    cout << "Trying out vectors\n";
    vector <Individual> prova;
    int dim = 4;
    int dim_seq = 3;
    prova.resize(dim);
    for (int i = 0; i< dim; i++) prova[i].setN(dim_seq);

    for (int i = 0; i<dim_seq; i++){ 
        prova[1].sequence.push_back(5);
        cout << prova[1].sequence[i] << " ";
    }
    
    cout << endl <<endl;
    prova.clear();

    for (int i = 0; i<dim_seq; i++){ 
        cout << prova[1].sequence[i] << " ";
    }


    cout << endl <<endl;
*/