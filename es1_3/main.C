#include "../lib.h"
#include "../random.h"

using namespace std;

int main(){
    ofstream out;
    out.open("es3.txt");

    Random rnd;
    rnd.Initialize(rnd);

    double d = 1; //distance among lines
    double s = d*0.7; //length of stick
    //creating environment. A 1D line is enough, since only horizontal direction is relevant in the intersection
    double xmin = 0, xmax = d*100;

    int M = 100000; //number of gen values
    int N = 100; //values per block
    int L = (int) M/N;
    int start = 1; //start from the block "start"

    for (int i = start; i < N; i++){ 
        int throws = (i+1)*L; //number of throws at this step
        int Nblocks  = i+1; //number of blocks at this step
        double * A = new double [Nblocks];
        for (int j = 0; j < Nblocks ; j++){
            double Nthr = 0, Nhit = 0;
            for (int k = 0; k < L; k++){
                Nthr++;
                double x = rnd.Rannyu(xmin,xmax); //center of stick
                double x_appo,y_appo;
                do{
                    x_appo = rnd.Rannyu();
                    y_appo = rnd.Rannyu();
                }while(pow(x_appo,2) + pow(y_appo,2) > 1.); //gen a random couple in a circumference of r=1
                double c = cos(atan(y_appo/x_appo)); //orientation in space 
                double pos = (x - int(x)); //position of the stick in the section between 2 lines
                double pro = s/2.*c; //projection of the stick on the x axis
                if (pro+pos > d or pro+pos < 0. or pos-pro < 0. or pos-pro > d){Nhit++;}
            }   
            double P = Nhit/Nthr;
            double pi = (2*s)/(P*d);
            A[j] = pi;
        }
        cout << "Blocks " << Nblocks << " Value: " << mean(Nblocks,A) <<endl;
        out << Nblocks << "," << mean(Nblocks,A) << "," << sqrt(variance_blocks(Nblocks,A)) << endl;
        delete []A;
    }

    out.close();
    return 0;
}
