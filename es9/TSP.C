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
        cout << "\n There's a problem with this element! First element is not first city" << endl;
        exit(0);
    }
    sort(v.begin(), v.end());
    if (v != ordered){
        cout << "\n There's a problem with this element!" << endl;
        for (int i = 0; i<v.size(); i++) cout << v[i]<<",";
        cout << endl;
        exit(0);
    }
}

double Genetic :: loss (vector <int> v){
    double sum = 0.;
    for (int i = 1; i < v.size(); i++) sum += pow(abs(cities[v[i]].Distance(cities[v[i-1]])),L);
    //add the loss between first and last
    sum += pow(abs(cities[v[0]].Distance(cities[v[v.size()-1]])),L); 
    return sum;
}

void Genetic :: order (vector <Individual> &a){
    int dim = a.size();
    for (int i = 0; i < dim; i++) 
        for (int j = i+1; j < dim; j++) 
            if (a[i].fitness > a[j].fitness) swap(a[i],a[j]);
}

void Genetic :: selection (vector <Individual> old_gen, Individual & new_ind){
    int new_index = (int) n_elems *(pow(rnd.Rannyu(), p));
    //cout << "\nIndex selected from previous generation is " << new_index << endl; 
    new_ind = old_gen[new_index];
    //cout << "Confronto le fitness: old = " << old_gen[new_index].fitness << " new = " << new_ind.fitness << endl;
}

void Genetic :: newGen (vector <Individual> previous_gen, vector <Individual> & next_gen){
    for (int i = 0; i < previous_gen.size(); i++) selection(previous_gen, next_gen[i]); 
}

void Genetic :: randomSwap (vector <Individual>& v){
    int a,b;
    for (int i = 1; i < v.size(); i++){ //leave the first element untouched
        if (rnd.Rannyu() < prob_mut){ 
            a = rnd.Rannyu(1,n_elems);
            do{b = rnd.Rannyu(1,n_elems);} while(a==b);
            swap (v[i].sequence[a],v[i].sequence[b]);
            check(v[i].sequence);
            v[i].fitness = loss(v[i].sequence);
        }
    }
    order(v);
}

void Genetic :: contSwap (vector <Individual>& v){
    //selecting elems in the middle to avoid going over dimension
  
    for (int i = 1; i < v.size(); i++){ //leave the first element untouched
        if (rnd.Rannyu() < prob_mut){
            int length = rnd.Rannyu(1,n_elems/2.);
            int which1 = rnd.Rannyu(1,n_elems- length), which2 = rnd.Rannyu(1,n_elems-length); 
            for (int k = 0; k < length; k++)  swap (v[i].sequence[which1 +k],v[i].sequence[which2 +k]);
            check(v[i].sequence);
            v[i].fitness = loss(v[i].sequence);
        }
    }
    order(v);
}

//this shifts a random element to the end of the sequence
void Genetic :: shift (vector <Individual>& v){
    for (int i = 1; i < v.size(); i++){ 
        if (rnd.Rannyu() < prob_mut){ //do this with 10% prob
            unsigned int which = (int) rnd.Rannyu(1,n_elems);
            v[i].sequence.push_back(v[i].sequence[which]);
            v[i].sequence.erase(v[i].sequence.begin()+which);
            check(v[i].sequence);
            v[i].fitness = loss(v[i].sequence);
        }
    }
    order(v);
}

void Genetic:: crossover(vector<Individual> &v){
    int j, start;
    for (int i = 1; i < v.size(); i++){ 
        if (rnd.Rannyu() < 0.5){
            do{j = rnd.Rannyu(1,v.size());} while(j==i);
            vector <int> seqi = v[i].sequence;
            vector <int> seqj = v[j].sequence;
            vector <int> newi = seqi;
            vector <int> newj = seqj;
            
            start = int(rnd.Rannyu(1,seqi.size())); //start of cross
            vector <int > indexi;
            vector <int > indexj;
            vector<int>::iterator it;
            for(int k = start; k < seqi.size(); k++){
                it = find (seqj.begin(), seqj.end(),seqi[k]);
                indexi.push_back(it-seqj.begin());
                it = find (seqi.begin(), seqi.end(),seqj[k]);
                indexj.push_back(it-seqi.begin());
            }
            sort(indexi.begin(),indexi.end());
            sort(indexj.begin(),indexj.end());
            //cout << "\nordered vector\n";
            //for(int k = 0; k < indexi.size();k++) cout << indexi[k] << ",";

            int t = 0;
            for(int k = start; k < seqi.size();k++){
                newi[k] = seqj[indexi[t]];
                newj[k] = seqi[indexj[t]];
                t++;
            }
            //cout << "\nold sequence is: \n";
            //for(int k = 0; k < seqi.size();k++) cout << v[i].sequence[k] << ",";
            v[i].sequence = newi;
            v[j].sequence = newj;        
            //cout << "\nnew sequence is: \n";            
            //for(int k = 0; k < seqi.size(); k++) cout << newi[k] << ",";
            //cout << endl;
            check(v[i].sequence);
            check(v[j].sequence);
            v[i].fitness = loss(v[i].sequence);
            v[j].fitness = loss(v[j].sequence);
        }
    }
    order(v);
}

void Genetic :: end (){
    rnd.SaveSeed();
}

