#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){
    ofstream out;
    out.open("es3.txt");

    Random rnd;
    rnd.Initialize(rnd);

    //creating environment. A 1D line is enough, since only horizontal direction is relevant in the intersection
    double d = 1;       //distance among lines
    double s = d*0.7;   //length of the stick
    double xmin = 0, xmax = d*100;

    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    int Nthr = 0, Nhit = 0;
    double A[N] = {};
    double x = 0, x_appo = 0, y_appo = 0, c = 0, pos = 0, pro = 0;
    double P = 0, pi = 0, var = 0;

    for (int i = 0; i < N; i++){ 
        int throws = (i+1)*L;   //number of throws at this step
        int Nblocks  = i+1;     //number of blocks at this step
        for (int j = 0; j < Nblocks ; j++){
            Nthr = 0;
            Nhit = 0;
            for (int k = 0; k < L; k++){
                Nthr++;
                x = rnd.Rannyu(xmin,xmax); //position of the center of the stick
                // generate a random couple of coor in a circumference of r=1, using hit or miss
                
                do{
                    x_appo = rnd.Rannyu();
                    y_appo = rnd.Rannyu();
                } while (pow(x_appo,2) + pow(y_appo,2) > 1.); 

                c = cos(atan(y_appo/x_appo)); //orientation in space 
                pos = (x - int(x));           //position of the stick in the section between 2 lines
                pro = s/2.*c;                 //projection of the stick on the x axis
                // now, count as hit every possible intersection case
                if (pro+pos > d or pro+pos < 0. or pos-pro < 0. or pos-pro > d){Nhit++;}
            }   
            P = (double) Nhit/Nthr;
            pi = (2*s)/(P*d);
            A[i] += pi;
        }
        A[i] /= (double) L;
        if (Nblocks != 1) var = variance_blocks(Nblocks,A);
        cout << "Block " << Nblocks << ". Value: " << A[i] <<endl;
        out << Nblocks << "," << mean(Nblocks, A) << "," << var << endl;
    }

    out.close();
    return 0;
}
