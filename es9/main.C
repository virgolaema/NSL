#include "TSP.h"

using namespace std;

int main(){
    int n = 32; //number of cities
    int n_gens = 32*32; //number of generations

    Genetic circ (n);
    circ.dispose(0); //dispose the n cities on a circumference
    circ.whichLoss(1); //1 for L1, 2 for L2

    ofstream pos ("circ.out");
    for (int i = 0; i < n; i++)
    pos << circ.cities[i].getX() << "," << circ.cities[i].getY() << endl;
    pos.close();

    cout << "Elements disposed on a circ \n";

    vector <Individual> first_gen; 
    first_gen.resize(n);

    cout << "Generated array of individuals, first gen \n";

    for (int i = 0; i < n; i++){
        first_gen[i].setN(n);
        circ.newIndividual(first_gen[i].sequence); //start creating random individuals
        circ.check(first_gen[i].sequence);
        first_gen[i].fitness = circ.loss(first_gen[i].sequence);
    }
    cout << endl << endl;

    // orderING BASING ON LOSS FUNC 1 or 2
    circ.order(first_gen);

    /*for (int i = 0; i < n; i++){
        cout << "\n The first created generation is:";
        cout <<"Loss: " << first_gen[i].fitness << endl;
        for (int h = 0; h < n; h++) cout <<  first_gen[i].sequence[h] << " ";
        cout << endl;
    }*/
    
    for (int i = 0; i < n; i++) circ.check(first_gen[i].sequence);

    vector <Individual> next_gen; 
    next_gen.resize(n);

    cout << "Created new vector for next gens\n";

    vector <Individual> previous_gen = first_gen;

    for (int i = 2; i < n_gens; i++){

        cout << "Creating the " << i << " generation...";
        circ.shift(previous_gen);
        circ.randomSwap(previous_gen);
        circ.order(previous_gen);
        circ.newGen(previous_gen, next_gen); //check included
        cout << "Done\n";
        circ.order(next_gen);
        cout << "Ordered\n";
        previous_gen = next_gen;
        cout << "Update "; 
    }

    for (int i = 0; i < n; i++){
        cout << "\n The last created generation is: \n";
        for (int h = 0; h < n; h++) cout << next_gen[i].sequence[h] << " ";
        cout << endl;
    }

    ofstream final("final.out");

    for (int h = 0; h < n; h++) 
        final << circ.cities[next_gen[0].sequence[h]-1].getX() << "," << circ.cities[next_gen[0].sequence[h]-1].getY() << endl;

    final.close();

    cout << circ.loss(next_gen[0].sequence);

    //now, get busy. How many generations? selection mutation


   Genetic sqre (n);
    sqre.dispose(1);


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