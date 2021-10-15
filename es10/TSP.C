#include "TSP.h"

using namespace std;


SimAnn :: SimAnn (int N):cities(N){
    n_elems = N;
    Random appo;
    rnd = appo;
    rnd.Initialize(rnd);
    //creating the ordered 0-31 vector, used in check()
    for (int i = 0; i < n_elems; i++)  ordered.push_back(i); 
}

//simple shuffle, every element gets swapped with another random one
// note that the first is always 1
void SimAnn:: shuffle(vector<int> &v){
    int j;
    for (int i = 0; i < v.size(); i++){
        j = 0;
        do{j = rnd.Rannyu(1,n_elems);} while(i == j); //avoid swapping with self
        swap(v[i], v[j]);
    }
}

void SimAnn :: dispose (int which){
    if (which == 0){ //circumference
        double angle = 0;
        for (int i = 0; i < n_elems; i++){
            //angle = 2.* (double)M_PI *(double) i / (double) n_cities; //this is to put them uniformly
            angle = rnd.Rannyu(0.,2*M_PI);  
            cities[i].setX(cos(angle));
            cities[i].setY(sin(angle));
        }
        cout << "Elements disposed on a circumference \n";

    }
    else if (which == 1)    { //inside a square
        for (int i = 0; i < n_elems; i++){
            cities[i].setX(rnd.Rannyu (-1.,1.));
            cities[i].setY(rnd.Rannyu (-1.,1.));
        }
        cout << "Elements disposed in a square \n";
    }
    else cout << "Choose a valid disposition: 0 for circumference, 1 for square" << endl;
}

//creating a new individual, a simple ordered vector from 1 to n_elems and then shuffle
void SimAnn:: newIndividual(vector<int> &v){
    for (int i = 0; i < n_elems; i++) v.push_back(i); //filled in order, from 1 to 32
    shuffle(v); //then shuffle to randomize it
    check(v);
}

// 1: check if the first element is the city #1
// 2: Puts the vector in order with order, then compare with an ordered one
void SimAnn :: check(vector <int> v){
   /* if (v[0] != 0){
        cout << "\n There's a problem with this element! First element is not first city" << endl;
        for (int i = 0; i<v.size(); i++) cout << v[i]<<",";
        exit(0);
    }*/
    sort(v.begin(), v.end());
    if (v != ordered){
        cout << "\n There's a problem with this element!" << endl;
        for (int i = 0; i<v.size(); i++) cout << v[i]<<",";
        cout << endl;
        exit(0);
    }
}

double SimAnn :: loss (vector <int> v){
    double sum = 0.;
    for (int i = 1; i < v.size(); i++) sum += pow(abs(cities[v[i]].Distance(cities[v[i-1]])),L);
    //add the loss between first and last
    sum += pow(abs(cities[v[0]].Distance(cities[v[v.size()-1]])),L); 
    return sum;
}

void SimAnn :: randomSwap (Individual & v){
    int a,b;
    if (rnd.Rannyu() < prob_mut){ 
        a = rnd.Rannyu(0,n_elems);
        do{b = rnd.Rannyu(0,n_elems);} while(a==b);
        swap (v.sequence[a],v.sequence[b]);
        check(v.sequence);
        v.fitness = loss(v.sequence);
    }
}

void SimAnn :: contSwap (Individual & v){
    if (rnd.Rannyu() < prob_mut){
        int length = rnd.Rannyu(0,n_elems/2.);
        int which1 = rnd.Rannyu(0,n_elems);
        int which2 = Pbc(which1+length); 
        for (int k = 0; k < length; k++)  swap (v.sequence[Pbc(which1 +k)],v.sequence[Pbc(which2 +k)]);
        check(v.sequence);
        v.fitness = loss(v.sequence);
    }
}

//this shifts a random element to the end of the sequence
void SimAnn :: shift (Individual & v){
    if (rnd.Rannyu() < prob_mut){ 
        unsigned int which = (int) rnd.Rannyu(0,n_elems);
        v.sequence.push_back(v.sequence[which]);
        v.sequence.erase(v.sequence.begin()+which);
        check(v.sequence);
        v.fitness = loss(v.sequence);
    }
}

int SimAnn :: Pbc (int a){ 
    if (a >= n_elems) a = 1 + a - n_elems;
    return a;
}

void SimAnn :: end (){rnd.SaveSeed();}

