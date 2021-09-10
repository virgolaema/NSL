#include "TSP.h"

using namespace std;

Genetic :: Genetic (int N){
    n_elems = N;
    Random appo;
    rnd = appo;
    rnd.Initialize(rnd);
    matrix = new int* [n_elems];
    for (int i = 0; i < n_elems; i++)
        matrix[i] = new int [n_elems];
    cities = new Position [n_elems];
}

void Genetic :: set1 (int row, int column){matrix[row][column] = 1;}

void Genetic :: set0 (int row, int column){matrix[row][column] = 0;}

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

void Genetic :: swap (int r1, int c1,int r2, int c2){
    int appo = matrix[r1][c1];
    matrix[r1][c1] = matrix[r2][c2];
    matrix[r2][c2] = appo;
}

void Genetic :: check(){
    //checking is there's only one value 1 in the columns or rows
    for (int i = 0; i < n_elems; i++){
        int sum_column = 0;
        int sum_row = 0;
        for (int j = 0; j< n_elems; j++){
            sum_column += matrix [j][i];
            sum_row += matrix [i][j];
        }
        if (sum_column != 1 or sum_row != 1){
            cout << "ALERT!!! Error in the matrix" << endl;
            return;
        }
    }
}

double Genetic :: loss1 (){
    double sum = 0.;
    for (int i = 0; i<n_elems; i++){
        int next_ = i+1;
        if (i = n_elems-1) next_ = 0;
        sum += abs(cities[i].Distance(cities[next_]));
    }
    return sum;
}

double Genetic :: loss2 (){
    double sum = 0.;
    for (int i = 0; i<n_elems; i++){
        int next_ = i+1;
        if (i = n_elems-1) next_ = 0;
        sum += pow(cities[i].Distance(cities[next_]),2);
    }
    return sum;
}

void Genetic :: end (){
    rnd.SaveSeed();
}