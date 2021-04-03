#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){
    double L = d*0.7; //length of stick
    double d = 1; //distance among lines
    double Nthr = 0, Nhit = 0;
    double P = Nhit/Nthr;

    //creating environment, a square 100x100
    double xmax = d*100, xmin = 0, ymin = 0; ymax = d*100;

    double x = randUnif(xmin,xmax); //center of stick
    double y = randUnif(ymin,ymax);
    double c = 



    double pi = (2*L)/(P*d);


    return 0;
}
