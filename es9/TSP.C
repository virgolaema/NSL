#include "TSP.h"

using namespace std;

void Individual :: setN (int n){
    N = n; 
    //sequence.resize(N);
    fitness = 0; //just to initialize
}

//////////////////////////////////////////////////////////////////////


Genetic :: Genetic (int N){
    n_elems = N;
    Random appo;
    rnd = appo;
    rnd.Initialize(rnd);
    cities = new Position [n_elems];
    //creating the ordered 1-32 vector, used in check()
    for (int i = 0; i < n_elems; i++)   ordered.push_back(i+1); 

}

//simple shuffle, every element gets swapped with another random one
// note that the first is always 1
void Genetic:: shuffle(vector<int> &v){
    for (int i = 1; i < v.size(); i++){
        int j = 0;
        do{j = (int) rnd.Rannyu(1,n_elems);} while(i == j); //avoid swapping with self
        swap(v[i], v[j]);
    }
}

void Genetic :: dispose (int which){
    if (which == 0){ //circumference
        double angle = 0;
        for (int i = 0; i < n_elems; i++){
            //angle = 2.* (double)M_PI *(double) i / (double) n_cities; //this is to put them uniformly
            angle = rnd.Rannyu(0.,2*M_PI);
            cities[i].setCoor(cos(angle), sin (angle),0.);
        }
    }
    else if (which == 1)    { //inside a square
        for (int i = 0; i < n_elems; i++)
            cities[i].setCoor(rnd.Rannyu (-1.,1.),rnd.Rannyu (-1.,1.),0.);
    }
    else cout << "Choose a valid disposition: 0 for circumference, 1 for square" << endl;
}

//creating a new individual, a simple ordered vector from 1 to n_elems and then shuffle
void Genetic:: newIndividual(vector<int> &v){
    for (int i = 0; i < n_elems; i++) v.push_back(i+1); //filled in order, from 1 to 32
    shuffle(v); //then shuffle to randomize it
}

// Checks if it's all good. 
// 1: check if the first element is the city #1
// 2: Puts the vector in order with order, then compare with an ordered one
void Genetic :: check(vector <int> v){
    if (v[0] != 1){
        cout << "\n There's a problem with this element!" << endl;
        exit(0);
    }
    sort(v.begin(), v.end());
    if (v != ordered){
        cout << "\n There's a problem with this element!" << endl;
        exit(0);
    }
}

double Genetic :: loss (vector <int> v){
    double sum = 0.;
    int next_ = 1;
    for (int i = 0; i<n_elems; i++){
        if (i = n_elems-1) next_ = 0;
        sum += abs(pow(cities[v[i]].Distance(cities[v[next_]]),L));
        next_++;
    }
    return sum;
}

void Genetic :: order (vector <Individual> &a){
    int dim = a.size();
    for (int i = 0; i < dim; i++) 
        for (int j = i+1; j < dim; j++) 
            if (a[i].fitness > a[j].fitness){
                Individual appo;
                appo = a[i];
                a[i] = a[j];
                a[j] = appo;
            }
}


void Genetic :: selection (vector <Individual> a, Individual & new_ind){
    int new_index = (int) n_elems *(1- pow(rnd.Rannyu(), 0.9));
    //cout << "\nIndex selected from previous generation is " << new_index << endl; 
    new_ind.sequence = a[new_index].sequence;
}

void Genetic :: newGen (vector <Individual> previous_gen, vector <Individual> & next_gen){
    next_gen.resize(n_elems);

    for (int i = 0; i < n_elems; i++){
        next_gen[i].setN(n_elems);
        //now,create the next generation from the previous
        selection(previous_gen, next_gen[i]); 
        check(next_gen[i].sequence);
        next_gen[i].fitness = loss(next_gen[i].sequence);
    }
}

void Genetic :: randomSwap (vector <Individual>& v){
    for (int i = 5; i < n_elems; i++) //leave the first 5 elements untouched, they are the best
        swap (v[i].sequence[rnd.Rannyu(1,n_elems)],v[i].sequence[rnd.Rannyu(1,n_elems)]);
    //might swap with itself, nevermind
}

//this shifts a random element to the end of the sequence
void Genetic :: shift (vector <Individual>& v){
    for (int i = 5; i < n_elems; i++){ //leave the first 5 elements untouched, they are the best
        int which = (int) rnd.Rannyu(1,n_elems);
        vector<int>::iterator it = v[i].sequence.begin() + which;
        v[i].sequence.push_back(v[i].sequence[which]);
        v[i].sequence.erase(it); 
    }
}

void Genetic :: end (){
    rnd.SaveSeed();
}

