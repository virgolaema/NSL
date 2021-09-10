#include "TSP.h"

using namespace std;

int main(){
    int n_cities = 32;

    Genetic circ (n_cities);
    circ.dispose(0);
    
    Genetic sqre (n_cities);
    sqre.dispose(1);





    /*ofstream disp_circ ("circ.out");
    ofstream disp_sqre ("sqre.out");

    for (int i = 0; i < n_cities; i++){
        disp_circ << circ.cities[i].getX() << "," << circ.cities[i].getY() << endl;
        disp_sqre << sqre.cities[i].getX() << "," << sqre.cities[i].getY() << endl;
    }

    disp_sqre.close();
    disp_circ.close();

    sqre.end();
    circ.end();*/

    return 0;
}